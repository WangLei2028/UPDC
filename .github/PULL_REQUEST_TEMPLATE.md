## Summary

Describe the change and why it is needed.

## Validation

- [ ] `python -m pytest -q`
- [ ] `updc doctor` checked when build/compiler behavior changed
- [ ] Documentation/config examples updated when interfaces changed
- [ ] Cache ABI/schema versions reviewed when persistent formats changed
- [ ] Frozen v17 C-core hashes are unchanged, **or** this PR explicitly proposes a new core freeze with scientific validation

## Scientific impact

State whether this change can alter reconstructed SCI/WHT values, WCS, PSF handling, weighting, background treatment, or iteration semantics.
