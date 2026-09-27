#!/usr/bin/env python3
"""Build the English UPDC manuals from Markdown with Pandoc + XeLaTeX.

Usage:
    python tools/build_manual_pdfs.py

Required developer tools: pandoc and xelatex.
Generated files are written next to their Markdown sources in docs/.
"""
from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"

MANUALS = [
    ("UPDC_COMPLETE_MANUAL_EN.md", "UPDC_COMPLETE_MANUAL_EN.pdf"),
    ("UPDC_ALGORITHM_MANUAL_EN.md", "UPDC_ALGORITHM_MANUAL_EN.pdf"),
    ("UPDC_MANUAL_EN.md", "UPDC_MANUAL_EN.pdf"),
]


def require(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"Required developer tool not found: {name}")
    return path


def build(md_name: str, pdf_name: str) -> None:
    pandoc = require("pandoc")
    require("xelatex")
    src = DOCS / md_name
    dst = DOCS / pdf_name
    cmd = [
        pandoc,
        str(src),
        "--from=gfm+tex_math_dollars",
        "--pdf-engine=xelatex",
        "--toc",
        "--toc-depth=2",
        "-V", "papersize=a4",
        "-V", "geometry:margin=22mm",
        "-V", "fontsize=10pt",
        "-V", "colorlinks=true",
        "-V", "linkcolor=blue",
        "-V", "urlcolor=blue",
        "-V", "monofont=DejaVu Sans Mono",
        "-o", str(dst),
    ]
    subprocess.run(cmd, check=True)
    if not dst.is_file() or dst.stat().st_size == 0:
        raise SystemExit(f"PDF generation failed: {dst}")
    print(dst)


def main() -> int:
    for md_name, pdf_name in MANUALS:
        build(md_name, pdf_name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
