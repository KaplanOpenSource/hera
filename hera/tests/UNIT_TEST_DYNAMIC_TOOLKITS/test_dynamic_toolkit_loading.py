import sys
import textwrap
from pathlib import Path
import pytest

# We test the new Class datatype handler
from hera.datalayer.datahandler import DataHandler_Class, getHandler


@pytest.fixture
def module_dir(tmp_path):
    """
    Create a temporary package 'mypkg' with module 'mymod.py' and class 'TestClass'.
    The handler is expected to add this directory to sys.path dynamically.
    We also clean sys.path after the test to avoid leaking state.
    """
    pkg = tmp_path / "mypkg"
    pkg.mkdir(parents=True, exist_ok=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "mymod.py").write_text(textwrap.dedent("""
        class TestClass:
            def __init__(self, a=0, b="ok"):
                self.a = a
                self.b = b
    """))
    # yield the directory that should be appended/inserted by the handler
    yield tmp_path
    # cleanup sys.path so other tests are not affected
    sp = str(tmp_path)
    while sp in sys.path:
        sys.path.remove(sp)


def test_instantiation_with_params(module_dir):
    """Handler returns an instance with constructor params applied."""
    obj = DataHandler_Class.getData(
        resource=str(module_dir),
        desc={
            "classpath": "mypkg.mymod.TestClass",
            "params": {"a": 7, "b": "wow"},
        },
    )
    assert obj.__class__.__name__ == "TestClass"
    assert (obj.a, obj.b) == (7, "wow")
    # the handler should have added the module_dir to sys.path
    assert str(module_dir) in sys.path


def test_via_getHandler(module_dir):
    """Same flow via the generic getHandler('Class')."""
    Handler = getHandler("Class")
    obj = Handler.getData(
        str(module_dir),
        desc={"classpath": "mypkg.mymod.TestClass", "params": {"a": 1, "b": "via-getHandler"}},
    )
    assert obj.__class__.__name__ == "TestClass"
    assert (obj.a, obj.b) == (1, "via-getHandler")


def test_return_class_not_instance(module_dir):
    """When instantiate=False, handler returns the class object (type) rather than an instance."""
    cls = DataHandler_Class.getData(
        resource=str(module_dir),
        desc={"classpath": "mypkg.mymod.TestClass", "instantiate": False},
    )
    assert isinstance(cls, type)
    assert cls.__name__ == "TestClass"


def test_errors_are_clear(module_dir):
    """Invalid inputs produce clear ValueError / ImportError."""
    # Missing classpath
    with pytest.raises(ValueError):
        DataHandler_Class.getData(resource=str(module_dir), desc={})

    # Bad module path
    with pytest.raises(ImportError):
        DataHandler_Class.getData(
            resource=str(module_dir),
            desc={"classpath": "mypkg.bad.DoesNotExist"}
        )

# ---------- New tests: DataHandler_Class ----------

from pathlib import Path
import textwrap
from hera.datalayer.datahandler import getHandler

def _write_pkg(root: Path):
    """Create a tiny temp package with a class to load dynamically."""
    (root / "mypkg").mkdir(parents=True, exist_ok=True)
    (root / "mypkg" / "__init__.py").write_text("")
    (root / "mypkg" / "mymod.py").write_text(textwrap.dedent("""
    class TestClass:
        def __init__(self, a=0, b="ok"):
            self.a = a
            self.b = b
        def __repr__(self):
            return f"<TestClass a={self.a} b={self.b}>"
    """))

def test_class_datatype_instance(tmp_path):
    """
    Verify DataHandler_Class returns an *instance* when instantiate=True (default),
    and passes constructor parameters correctly.
    """
    _write_pkg(tmp_path)
    resource = str(tmp_path)  # folder to be appended to sys.path by the handler
    desc = {
        "classpath": "mypkg.mymod.TestClass",
        "parameters": {"a": 7, "b": "wow"},
        # "instantiate": True  # default
    }
    handler = getHandler("Class")
    obj = handler.getData(resource, desc)
    assert obj.__class__.__name__ == "TestClass"
    assert (obj.a, obj.b) == (7, "wow")

def test_class_datatype_type_only(tmp_path):
    """
    Verify DataHandler_Class returns the *class* when instantiate=False.
    """
    _write_pkg(tmp_path)
    resource = str(tmp_path)
    desc = {
        "classpath": "mypkg.mymod.TestClass",
        "instantiate": False,
    }
    handler = getHandler("Class")
    cls = handler.getData(resource, desc)
    assert isinstance(cls, type)
    inst = cls(a=1, b="via-type")
    assert (inst.a, inst.b) == (1, "via-type")
