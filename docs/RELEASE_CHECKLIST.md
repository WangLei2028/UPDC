# Release checklist

- [ ] `updc doctor` reports all frozen-core files `OK`.
- [ ] `pytest -q` passes on supported Python versions.
- [ ] C core builds/self-tests on Linux and macOS with FFTW3f.
- [ ] Wheel builds and installs into a clean environment.
- [ ] `updc init`, `updc validate`, and `updc adapters` work from the installed wheel.
- [ ] Schema version unchanged unless an intentional migration exists.
- [ ] Cache ABI versions reviewed; incompatible format changes increment ABI.
- [ ] v17 reference benchmark compared when performance-relevant orchestration changes are made.
- [ ] Documentation and changelog updated.
- [ ] MIT `LICENSE`, `NOTICE`, `CITATION.cff`, and `CITATION.md` reviewed; both required UPDC papers are present.
