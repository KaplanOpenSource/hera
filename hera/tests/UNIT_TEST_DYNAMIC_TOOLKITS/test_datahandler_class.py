# Tests for DataHandler_Class merge rule (desc.parameters overrides kwargs)
import textwrap
from hera.datalayer.datahandler import DataHandler_Class

def _write_dummy_package(tmp_path):
    """
    Create a minimal package structure on disk:
      tmp_path/
        mypkg/
          __init__.py
          mycls.py   (defines class Echo)
    """
    pkg = tmp_path / "mypkg"
    pkg.mkdir(parents=True, exist_ok=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "mycls.py").write_text(textwrap.dedent("""
    class Echo:
        def __init__(self, a=0, b=0, **kwargs):
            self.a = a
            self.b = b
            self.kw = dict(kwargs)
    """))
    return tmp_path  # return the folder to add to sys.path via `resource`

def test_desc_parameters_override_kwargs(tmp_path):
    root = _write_dummy_package(tmp_path)
    # desc.parameters should override duplicates from **kwargs (Option B)
    desc = {"classpath": "mypkg.mycls.Echo", "parameters": {"a": 1}}
    obj = DataHandler_Class.getData(resource=str(root), desc=desc, a=9, b=2, extra="X")
    assert obj.a == 1            # desc wins over kwargs
    assert obj.b == 2            # came from kwargs
    assert obj.kw.get("extra") == "X"

def test_instantiate_false_returns_class(tmp_path):
    root = _write_dummy_package(tmp_path)
    cls = DataHandler_Class.getData(resource=str(root),
                                    desc={"classpath": "mypkg.mycls.Echo", "instantiate": False})
    assert isinstance(cls, type)
    inst = cls(a=3)
    assert inst.a == 3

def test_invalid_classpath_raises(tmp_path):
    root = _write_dummy_package(tmp_path)
    try:
        DataHandler_Class.getData(resource=str(root),
                                  desc={"classpath": "bad.missing.Class"})
    except Exception as e:
        # We accept either ImportError or ValueError from the handler
        assert isinstance(e, (ImportError, ValueError))
    else:
        assert False, "Expected an exception for invalid classpath"
