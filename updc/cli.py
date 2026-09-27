"""Stable UPDC command line interface."""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

from . import __version__
from .adapters import AdapterContext, adapter_descriptions, get_adapter
from .cache_abi import CACHE_ABI_SUMMARY
from .config import engine_argv, load_config, validate_config, validate_paths
from .frozen_core import assert_frozen_core, verify_frozen_core
from .build import compiler_identity, discover_fftw, native_arch_flag, probe_c11_atomics, sanitize_c_environment


def _template_path() -> Path:
    return Path(__file__).resolve().parent / "templates" / "updc.toml"


def validate(path: Path, overrides: list[str] | None = None, quiet: bool = False) -> dict:
    cfg = load_config(path, overrides or [])
    adapter = get_adapter(cfg["instrument"]["adapter"])
    cfg = adapter.apply_defaults(cfg)
    validate_config(cfg)
    errors = adapter.validate(cfg)
    if errors:
        raise ValueError("; ".join(errors))
    validate_paths(cfg)
    assert_frozen_core()
    if not quiet:
        print(f"Configuration valid: {Path(path).expanduser().resolve()}")
        print(f"Schema: v{cfg['updc']['schema_version']}  adapter: {adapter.name}")
        print("Frozen v17 C core integrity: OK")
    return cfg


def run(path: Path, overrides: list[str] | None = None) -> None:
    cfg = validate(path, overrides, quiet=True)
    adapter = get_adapter(cfg["instrument"]["adapter"])
    with tempfile.TemporaryDirectory(prefix="updc-adapter-") as td:
        context = AdapterContext(config_path=Path(path).expanduser().resolve(), staging_dir=Path(td))
        prepared = adapter.prepare(cfg, context)
        cfg2 = adapter.apply_defaults(prepared.config)
        validate_config(cfg2)
        errors = adapter.validate(cfg2)
        if errors:
            raise ValueError("; ".join(errors))
        validate_paths(cfg2)
        print(f"UPDC {__version__} | schema v{cfg2['updc']['schema_version']} | adapter {adapter.name}", flush=True)
        from . import _engine
        _engine.main(engine_argv(cfg2))


def init_config(path: Path, force: bool = False) -> None:
    target = path.expanduser().resolve()
    if target.exists() and not force:
        raise FileExistsError(f"Refusing to overwrite {target}; use --force")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(_template_path(), target)
    print(f"Wrote {target}")


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="updc", description="UPDC signed-x astronomical image reconstruction")
    p.add_argument("--version", action="version", version=f"UPDC {__version__}")
    sub = p.add_subparsers(dest="command")
    runp = sub.add_parser("run", help="run a reconstruction")
    runp.add_argument("config", type=Path); runp.add_argument("--set", action="append", default=[], metavar="SECTION.KEY=VALUE")
    valp = sub.add_parser("validate", help="validate configuration and frozen-core integrity")
    valp.add_argument("config", type=Path); valp.add_argument("--set", action="append", default=[], metavar="SECTION.KEY=VALUE")
    initp = sub.add_parser("init", help="write a schema-v1 configuration template")
    initp.add_argument("path", nargs="?", type=Path, default=Path("updc.toml")); initp.add_argument("--force", action="store_true")
    sub.add_parser("adapters", help="list installed instrument adapters")
    docp = sub.add_parser("doctor", help="show installation/core/cache/build diagnostics")
    docp.add_argument("--cc", default="cc", help="C compiler to probe (default: cc)")
    return p


def main(argv=None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    # Compact primary workflow: `updc config.toml [--set ...]`.
    if argv and argv[0] not in {"run", "validate", "init", "adapters", "doctor", "-h", "--help", "--version"}:
        argv.insert(0, "run")
    ns = _parser().parse_args(argv)
    try:
        if ns.command == "run": run(ns.config, ns.set)
        elif ns.command == "validate": validate(ns.config, ns.set)
        elif ns.command == "init": init_config(ns.path, ns.force)
        elif ns.command == "adapters":
            for name, description, source in adapter_descriptions(): print(f"{name:16} {source:7} {description}")
        elif ns.command == "doctor":
            print(f"UPDC {__version__}")
            print("Frozen core:")
            for name, ok in verify_frozen_core().items(): print(f"  {name:20} {'OK' if ok else 'MISMATCH'}")
            print("Cache ABI:")
            for name, version in CACHE_ABI_SUMMARY.items(): print(f"  {name:20} v{version}")
            print("Build environment:")
            benv, removed = sanitize_c_environment()
            ident = compiler_identity(ns.cc, benv)
            print(f"  compiler             {ident['path']}")
            print(f"  compiler version     {ident['version']}")
            if removed:
                print("  C++ header shadow    detected and filtered for C builds")
                for key, paths in sorted(removed.items()):
                    print(f"    {key}: {os.pathsep.join(paths)}")
            else:
                print("  C++ header shadow    none detected")
            ok, _ = probe_c11_atomics(ident['path'], benv)
            print(f"  C11 atomics          {'OK' if ok else 'FAILED'}")
            native = native_arch_flag(ident['path'], benv)
            print(f"  native CPU flag      {native or 'unavailable'}")
            try:
                fftw = discover_fftw(ident['path'], benv)
                print(f"  FFTW3f               OK ({fftw['version']}, {fftw['source']})")
            except RuntimeError as exc:
                print(f"  FFTW3f               FAILED ({str(exc).splitlines()[0]})")
            print("Adapters:")
            for name, description, source in adapter_descriptions(): print(f"  {name:16} {source:7} {description}")
        else:
            _parser().print_help()
    except (ValueError, KeyError, FileExistsError, RuntimeError, OSError) as exc:
        raise SystemExit(f"updc: error: {exc}") from exc
