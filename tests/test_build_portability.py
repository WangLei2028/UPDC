from pathlib import Path

from updc.build import (
    build_signature, native_arch_flag, probe_c11_atomics, sanitize_c_environment,
)


def test_sanitize_filters_cpp_header_shadowing():
    env = {
        "PATH": "/usr/bin",
        "C_INCLUDE_PATH": "/opt/gcc/include/c++/13:/opt/site/include",
        "CPLUS_INCLUDE_PATH": "/opt/gcc/include/c++/13",
        "CPATH": "/a/include:/b/include/c++/13/x86_64:/c/include",
    }
    clean, removed = sanitize_c_environment(env)
    assert clean["C_INCLUDE_PATH"] == "/opt/site/include"
    assert clean["CPATH"] == "/a/include:/c/include"
    assert "CPLUS_INCLUDE_PATH" not in clean
    assert removed["C_INCLUDE_PATH"] == ["/opt/gcc/include/c++/13"]


def test_build_signature_changes_with_toolchain():
    a = build_signature("abc", {"path": "/gcc", "version": "gcc 13", "requested": "gcc"}, ["-O3"])
    b = build_signature("abc", {"path": "/gcc", "version": "gcc 14", "requested": "gcc"}, ["-O3"])
    c = build_signature("abc", {"path": "/gcc", "version": "gcc 13", "requested": "gcc"}, ["-O2"])
    assert a["signature_sha256"] != b["signature_sha256"]
    assert a["signature_sha256"] != c["signature_sha256"]


def test_host_c_compiler_supports_c11_atomics():
    clean, _ = sanitize_c_environment()
    ok, output = probe_c11_atomics("cc", clean)
    assert ok, output
    # Native flag is optional but the probe itself must not raise.
    native_arch_flag("cc", clean)
