"""Portable C build helpers for UPDC's frozen numerical core.

This module is packaging/build infrastructure only.  It intentionally does not
change the frozen numerical C sources.  Its goals are to make runtime builds
reproducible across macOS/Linux/HPC environments and to prevent stale compiler
binaries from being reused after the toolchain changes.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable

_CPP_INCLUDE_RE = re.compile(r"(?:^|[/\\])include[/\\]c\+\+(?:[/\\]|$)")


def _split_paths(value: str | None) -> list[str]:
    if not value:
        return []
    return [p for p in value.split(os.pathsep) if p]


def sanitize_c_environment(env: dict[str, str] | None = None) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Return an environment safe for compiling C sources.

    Some HPC compiler modules put libstdc++ directories such as
    ``.../include/c++/13`` into ``C_INCLUDE_PATH``.  That can shadow GCC's C11
    ``<stdatomic.h>`` with the C++ compatibility header when compiling C.  C++
    include paths are irrelevant to UPDC's C build, so they are filtered out.

    ``CPLUS_INCLUDE_PATH`` is removed entirely for the C compiler.  Other C
    include paths (for site software, FFTW, etc.) are preserved.
    """
    out = dict(os.environ if env is None else env)
    removed: dict[str, list[str]] = {}

    cplus = _split_paths(out.pop("CPLUS_INCLUDE_PATH", None))
    if cplus:
        removed["CPLUS_INCLUDE_PATH"] = cplus

    for key in ("CPATH", "C_INCLUDE_PATH"):
        parts = _split_paths(out.get(key))
        if not parts:
            continue
        keep: list[str] = []
        drop: list[str] = []
        for part in parts:
            normalized = part.replace("\\", "/")
            if _CPP_INCLUDE_RE.search(normalized):
                drop.append(part)
            else:
                keep.append(part)
        if drop:
            removed[key] = drop
            if keep:
                out[key] = os.pathsep.join(keep)
            else:
                out.pop(key, None)
    return out, removed


def compiler_identity(cc: str, env: dict[str, str] | None = None) -> dict[str, str]:
    env = os.environ if env is None else env
    resolved = shutil.which(cc, path=env.get("PATH")) or cc
    try:
        real = str(Path(resolved).expanduser().resolve())
    except OSError:
        real = str(resolved)
    try:
        text = subprocess.check_output([resolved, "--version"], text=True, stderr=subprocess.STDOUT, env=env)
        version = text.splitlines()[0].strip() if text else "unknown"
    except (OSError, subprocess.CalledProcessError):
        version = "unknown"
    return {"requested": str(cc), "path": real, "version": version}


def _compile_probe(cc: str, source: str, flags: Iterable[str] = (), libs: Iterable[str] = (),
                   env: dict[str, str] | None = None) -> tuple[bool, str]:
    env = os.environ if env is None else env
    with tempfile.TemporaryDirectory(prefix="updc-cc-probe-") as td:
        src = Path(td) / "probe.c"
        exe = Path(td) / "probe"
        src.write_text(source)
        cmd = [cc, *flags, str(src), "-o", str(exe), *libs]
        try:
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, env=env)
        except OSError as exc:
            return False, str(exc)
        return proc.returncode == 0, proc.stdout


def probe_c11_atomics(cc: str, env: dict[str, str] | None = None) -> tuple[bool, str]:
    source = r'''#include <stdint.h>
#include <stdatomic.h>
static atomic_uint_fast64_t x;
int main(void) {
    atomic_init(&x, 0);
    atomic_fetch_add_explicit(&x, 1, memory_order_relaxed);
    return atomic_load_explicit(&x, memory_order_relaxed) == 1 ? 0 : 1;
}
'''
    return _compile_probe(cc, source, flags=["-std=c11"], env=env)


def native_arch_flag(cc: str, env: dict[str, str] | None = None) -> str | None:
    """Return a supported native-CPU flag, or None.

    The flag is probed independently so unrelated source/link errors are never
    mislabeled as a rejection of ``-march=native``/``-mcpu=native``.
    """
    if sys.platform == "darwin" and platform.machine().lower() in {"arm64", "aarch64"}:
        candidate = "-mcpu=native"
    else:
        candidate = "-march=native"
    ok, _ = _compile_probe(cc, "int main(void){return 0;}\n",
                           flags=[candidate, "-std=c11"], env=env)
    return candidate if ok else None


def _pkg_config(args: list[str], env: dict[str, str]) -> str:
    return subprocess.check_output(["pkg-config", *args], text=True,
                                   stderr=subprocess.STDOUT, env=env).strip()


def _manual_fftw_candidates(env: dict[str, str]) -> tuple[list[Path], list[Path]]:
    incs: list[Path] = []
    libs: list[Path] = []
    if env.get("UPDC_FFTW_INCLUDE"):
        incs.append(Path(env["UPDC_FFTW_INCLUDE"]).expanduser())
    if env.get("UPDC_FFTW_LIB"):
        libs.append(Path(env["UPDC_FFTW_LIB"]).expanduser())
    for prefix in filter(None, [env.get("CONDA_PREFIX"), sys.prefix]):
        p = Path(prefix)
        incs.append(p / "include")
        libs.append(p / "lib")
    incs.extend([Path("/usr/include"), Path("/usr/local/include")])
    machine = platform.machine().lower()
    multiarch = {
        "x86_64": "x86_64-linux-gnu",
        "amd64": "x86_64-linux-gnu",
        "aarch64": "aarch64-linux-gnu",
        "arm64": "aarch64-linux-gnu",
    }.get(machine)
    libs.extend([Path("/usr/lib64"), Path("/usr/local/lib64"), Path("/usr/lib"), Path("/usr/local/lib")])
    if multiarch:
        libs.insert(0, Path("/usr/lib") / multiarch)
    # Stable de-duplication.
    incs = list(dict.fromkeys(incs))
    libs = list(dict.fromkeys(libs))
    return incs, libs


def discover_fftw(cc: str, env: dict[str, str] | None = None) -> dict:
    """Discover single-precision FFTW + thread wrapper.

    Prefer pkg-config.  If metadata is missing, fall back to common system,
    Conda/micromamba, and explicit UPDC_FFTW_{INCLUDE,LIB} locations, then
    verify the result by compiling and linking a small FFTW-thread program.
    """
    env = os.environ if env is None else env
    try:
        cflags = _pkg_config(["--cflags", "fftw3f"], env).split()
        fftw_libs = _pkg_config(["--libs", "fftw3f"], env).split()
        version = _pkg_config(["--modversion", "fftw3f"], env) or "unknown"
        thread_ok = subprocess.run(["pkg-config", "--exists", "fftw3f_threads"],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                   env=env).returncode == 0
        if thread_ok:
            cflags += _pkg_config(["--cflags", "fftw3f_threads"], env).split()
            thread_libs = _pkg_config(["--libs", "fftw3f_threads"], env).split()
        else:
            thread_libs = ["-lfftw3f_threads"]
        libs = thread_libs + fftw_libs + ["-lm"]
        source = "pkg-config"
    except (FileNotFoundError, subprocess.CalledProcessError):
        incs, libdirs = _manual_fftw_candidates(env)
        inc = next((p for p in incs if (p / "fftw3.h").is_file()), None)
        if inc is None:
            raise RuntimeError(
                "FFTW3 single-precision development headers were not found. Install FFTW3 development files "
                "or set UPDC_FFTW_INCLUDE=/path/include and UPDC_FFTW_LIB=/path/lib."
            )
        # Prefer a directory that actually contains both base and thread libs;
        # otherwise rely on the system linker search path.
        def has_lib(d: Path, stem: str) -> bool:
            return any(d.glob(stem + ".so*")) or any(d.glob(stem + ".dylib")) or any(d.glob(stem + ".a"))
        libdir = next((p for p in libdirs if has_lib(p, "libfftw3f") and has_lib(p, "libfftw3f_threads")), None)
        cflags = [f"-I{inc}"]
        libs = ([f"-L{libdir}"] if libdir else []) + ["-lfftw3f_threads", "-lfftw3f", "-lm"]
        version = "unknown"
        source = "filesystem-fallback"

    probe = r'''#include <fftw3.h>
int main(void) {
    if (!fftwf_init_threads()) return 2;
    fftwf_plan_with_nthreads(1);
    fftwf_cleanup_threads();
    return 0;
}
'''
    ok, output = _compile_probe(cc, probe, flags=["-std=c11", *cflags], libs=libs, env=env)
    if not ok:
        raise RuntimeError(
            "FFTW3f was discovered but a compile/link probe failed. Check fftw3.h, libfftw3f, "
            "libfftw3f_threads, compiler/library ABI, and PKG_CONFIG_PATH/LD_LIBRARY_PATH.\n" + output[-4000:]
        )
    return {"cflags": cflags, "libs": libs, "version": version, "source": source}


def build_signature(source_digest: str, compiler: dict[str, str], flags: list[str],
                    libraries: list[str] | None = None, fftw: dict | None = None) -> dict:
    payload = {
        "schema": 1,
        "source_sha256": source_digest,
        "compiler": compiler,
        "flags": list(flags),
        "libraries": list(libraries or []),
        "fftw": ({"version": fftw.get("version"), "source": fftw.get("source"),
                  "cflags": fftw.get("cflags"), "libs": fftw.get("libs")} if fftw else None),
        "platform": {"system": platform.system(), "machine": platform.machine()},
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["signature_sha256"] = hashlib.sha256(blob).hexdigest()
    return payload


def signature_matches(path: Path, expected: dict) -> bool:
    try:
        saved = json.loads(path.read_text())
    except (OSError, ValueError, TypeError):
        return False
    return saved.get("signature_sha256") == expected.get("signature_sha256")


def write_signature(path: Path, signature: dict) -> None:
    path.write_text(json.dumps(signature, indent=2, sort_keys=True) + "\n")
