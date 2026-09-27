from updc.adapters import InstrumentAdapter, adapter_descriptions, get_adapter


def test_builtin_adapters_are_registered():
    assert isinstance(get_adapter("jwst"), InstrumentAdapter)
    assert isinstance(get_adapter("fits"), InstrumentAdapter)
    names = {row[0] for row in adapter_descriptions()}
    assert {"jwst", "fits"} <= names
