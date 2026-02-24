import pytest

def _find_getdata_callable():
    """
    Hera has a static getData(resource, desc, **kwargs) implementation in:
    hera/datalayer/datahandler.py
    We import that module and return the function safely.
    """
    from hera.datalayer import datahandler as dh
    # in your codebase getData is a @staticmethod inside a class in that module;
    # easiest: grab the function by name from the module-level class that defines it.
    # We'll scan for any attribute that has a getData callable with signature (resource, desc, **kwargs).
    for name in dir(dh):
        obj = getattr(dh, name)
        if hasattr(obj, "getData") and callable(getattr(obj, "getData")):
            return getattr(obj, "getData")
    raise RuntimeError("Could not locate a getData(...) callable in hera.datalayer.datahandler")

@pytest.mark.unit
def test_dynamic_class_loading_via_getdata(tmp_path):
    """
    Unit test:
    - create a tiny python package on the fly
    - dynamically load a class using Hera getData(resource=..., desc={'classpath': ...})
    No DB / no CLI.
    """
    pkg = tmp_path / "mypkg"
    pkg.mkdir()
    (pkg / "__init__.py").write_text("")
    (pkg / "impl.py").write_text(
        "class MyDynamicClass:\n"
        "    def __init__(self):\n"
        "        self.ok = True\n"
    )

    getData = _find_getdata_callable()

    obj = getData(
        resource=str(pkg),
        desc={"classpath": "mypkg.impl.MyDynamicClass"},
        # לא להעביר dataFormat פה!
    )

    assert obj.ok is True
