# GitHub release checklist

0. Confirm the repository-root `LICENSE` is the standard MIT License, and `CITATION.cff`, `CITATION.md`, and `NOTICE` contain both required UPDC paper citations.
1. Confirm frozen core integrity: `updc doctor` and `pytest -q`.
2. Run a representative warm-cache scientific reconstruction and archive timing/metrics.
3. Build: `python -m build`.
4. Validate distributions: `python -m twine check dist/*`.
5. Verify wheel contents contain `updc/c_core`, JSON schema, TOML template, cache/core manifests, and docs in the source distribution.
6. Generate SHA-256 checksums for source archives, wheel, manual PDF, and release metadata as appropriate.
7. Tag the release (`v1.1.0`) only after the GitHub Actions matrix passes.
8. Upload release assets: wheel, source archive, complete English PDF manual, algorithm PDF manual, and SHA-256 file.
9. Do not modify frozen v17 C sources in packaging-only releases. A numerical-core change requires a new freeze version and scientific validation.
