# -*- coding: utf-8 -*-
import os
import sys
import pathlib
import textwrap
import importlib
import pytest

from hera.datalayer.datahandler import DataHandler_Class, getHandler
from hera.toolkit import ToolkitHome
from hera.datalayer import Project

# ---------- Helpers ----------

def _write_pkg_structure(base: pathlib.Path, as_package: bool = True):
    """
    Create a minimal module structure under 'base':
    - If as_package=True:
        base/
          mypkg/
            __init__.py
            mymod.py  -> contains TestClass
      classpath: 'mypkg.mymod.TestClass'
      resource: base/'mypkg'  (package dir)
    - If as_package=False:
        base/
          mypkg/
            mymod.py  -> contains TestClass
      classpath: 'mypkg.mymod.TestClass'
      resource: base   (top dir, not a package)
    """
    base.mkdir(parents=True, exist_ok=True)
    pkg_dir = base / "mypkg"
    pkg_dir.mkdir(parents=True, exist_ok=True)
    if as_package:
        (pkg_dir / "__init__.py").write_text("", encoding="utf-8")

    (pkg_dir / "mymod.py").write_text(textwrap.dedent("""
    class TestClass:
        def __init__(self, a=0, b="ok", **kw):
            self.a = a
            self.b = b
            self.kw = kw
        def whoami(self):
            return f"TestClass(a={self.a}, b={self.b})"
    """).strip(), encoding="utf-8")

    classpath = "mypkg.mymod.TestClass"
    resource = str(pkg_dir if as_package else base)
    return resource, classpath


def _make_demo_toolkit(base: pathlib.Path):
    """
    Create a module with DemoToolkit class for registerToolkit flow.
    Returns (resource_dir, classpath).
    """
    pkg = base / "demo_pkg"
    pkg.mkdir(parents=True, exist_ok=True)
    (pkg / "__init__.py").write_text("", encoding="utf-8")
    (pkg / "demo.py").write_text(textwrap.dedent("""
    class DemoToolkit:
        def __init__(self, projectName=None, filesDirectory=None, alpha=1, **kw):
            self.projectName = projectName
            self.filesDirectory = filesDirectory
            self.alpha = alpha
        def ping(self):
            return f"OK:{self.projectName}:{self.alpha}"
    """).strip(), encoding="utf-8")
    return str(pkg), "demo_pkg.demo.DemoToolkit"


# ---------- Tests ----------

def test_sys_path_package_dir(tmp_path):
    """
    When 'resource' points to a *package* directory (has __init__.py),
    the handler should add *both* the package dir and its parent to sys.path,
    so imports like 'mypkg.mymod' resolve correctly.
    """
    resource, classpath = _write_pkg_structure(tmp_path, as_package=True)
    obj = DataHandler_Class.getData(resource=resource,
                                    desc={"classpath": classpath, "parameters": {"a": 1, "b": "X"}})
    assert obj.whoami() == "TestClass(a=1, b=X)"


def test_sys_path_non_package_dir(tmp_path):
    """
    When 'resource' is a top directory (not a package),
    adding just that directory should be sufficient for 'mypkg.mymod' to resolve.
    """
    resource, classpath = _write_pkg_structure(tmp_path, as_package=False)
    obj = DataHandler_Class.getData(resource=resource,
                                    desc={"classpath": classpath, "parameters": {"a": 2, "b": "Y"}})
    assert obj.whoami() == "TestClass(a=2, b=Y)"


def test_instantiate_false_returns_class(tmp_path):
    """
    With instantiate=False the handler returns the class (type) rather than an instance.
    """
    resource, classpath = _write_pkg_structure(tmp_path, as_package=True)
    cls = DataHandler_Class.getData(resource=resource, desc={"classpath": classpath, "instantiate": False})
    assert isinstance(cls, type)
    inst = cls(a=3, b="Z")
    assert inst.whoami() == "TestClass(a=3, b=Z)"


def test_bad_classpath_raises_clean_error(tmp_path):
    """
    A non-existing classpath should raise ImportError with a clear message.
    """
    resource, _ = _write_pkg_structure(tmp_path, as_package=True)
    with pytest.raises(ImportError):
        DataHandler_Class.getData(resource=resource, desc={"classpath": "mypkg.missing.Nope"})


def test_registerToolkit_and_getToolkit_fallback(tmp_path, monkeypatch):
    """
    End-to-end:
    - create DemoToolkit on disk
    - register it as a ToolkitDataSource
    - fetch it back via DataHandler_Class and via ToolkitHome.getToolkit
    """
    import sys
    import importlib
    from hera.toolkit import ToolkitHome
    from hera.datalayer.datahandler import DataHandler_Class

    project = "UnitTestProject"
    th = ToolkitHome()

    # Create a demo toolkit module on disk
    resource_dir, classpath = _make_demo_toolkit(tmp_path)

    # IMPORTANT: for `import demo_pkg...` we must add the *parent* of the package dir
    parent_dir = os.path.dirname(resource_dir)
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)

    # Import the class object to pass into registerToolkit
    mod_name, _, cls_name = classpath.rpartition(".")
    mod = importlib.import_module(mod_name)
    DemoToolkit = getattr(mod, cls_name)

    # Register via ToolkitHome.registerToolkit
    doc = th.registerToolkit(
        toolkitclass=DemoToolkit,
        datasource_name="DemoToolkit_DS",
        params={"alpha": 7},
        version=(0, 0, 1),
        overwrite=True,
        projectName=project,
    )

    # Basic assertions on saved document
    assert isinstance(doc.resource, str) and os.path.isdir(doc.resource)
    assert doc.desc.get("classpath") == classpath

    # Instantiate via DataHandler_Class
    obj = DataHandler_Class.getData(resource=doc.resource, desc=doc.desc)
    assert obj.__class__.__name__ == "DemoToolkit"
    assert getattr(obj, "alpha", None) == 7

    # Instantiate via ToolkitHome.getToolkit (DB fallback)
    obj2 = th.getToolkit("DemoToolkit_DS", projectName=project)
    assert obj2.__class__.__name__ == "DemoToolkit"


def test_toolkit_table_contains_dynamic(tmp_path):
    """
    After registering a dynamic toolkit, getToolkitTable(project) should return
    a row for it, include 'source' column, and have the expected classpath.
    """
    import sys
    import importlib
    from hera.toolkit import ToolkitHome

    project = "UnitTestProject"
    th = ToolkitHome()

    # Create + import a class from a temp package
    resource_dir, classpath = _make_demo_toolkit(tmp_path)

    # Add PARENT of the package dir to sys.path so `import demo_pkg...` works
    parent_dir = os.path.dirname(resource_dir)
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)

    mod_name, _, cls_name = classpath.rpartition(".")
    DemoToolkit = getattr(importlib.import_module(mod_name), cls_name)

    # Register
    th.registerToolkit(
        toolkitclass=DemoToolkit,
        datasource_name="DemoToolkit_DS",
        params={"alpha": 3},
        version=(0, 0, 1),
        overwrite=True,
        projectName=project,
    )

    # Table should include our dynamic toolkit
    df = th.getToolkitTable(project)
    assert "toolkit" in df.columns
    assert "cls" in df.columns
    assert "source" in df.columns  # column exists
    assert (df["toolkit"] == "DemoToolkit_DS").any()

    row = df.loc[df["toolkit"] == "DemoToolkit_DS"].iloc[0]
    assert isinstance(row["cls"], str) and row["cls"].endswith(".DemoToolkit")


# ---------------- Optional smoke tests (skipped unless env is configured) ----------------

@pytest.mark.skipif(
    not (pathlib.Path(os.path.expanduser("~/hera-ims")).exists() and
         pathlib.Path(os.path.expanduser("~/hera-ims/token.json")).exists()),
    reason="IMS repo/token not available on this machine"
)
def test_ims_smoke_loads_from_repository():
    project = "UnitTestProject"
    p = Project(projectName=project)
    docs = p.getMeasurementsDocuments(type="ToolkitDataSource", datasourceName="IMS")
    assert docs, "IMS datasource not registered"
    ims = DataHandler_Class.getData(resource=docs[0].resource, desc=docs[0].desc)
    assert hasattr(ims, "download") or hasattr(ims, "update")


@pytest.mark.skipif(
    not pathlib.Path(os.path.expanduser("~/hera-jerusalem-main")).exists(),
    reason="Jerusalem repo not available on this machine"
)
def test_jerusalem_smoke_loads_from_repository():
    project = "UnitTestProject"
    p = Project(projectName=project)
    docs = p.getMeasurementsDocuments(type="ToolkitDataSource", datasourceName="Jerusalem2018_EXT")
    assert docs, "Jerusalem2018_EXT datasource not registered"
    obj = DataHandler_Class.getData(resource=docs[0].resource, desc=docs[0].desc)
    # just verify instance creation
    assert obj.__class__.__name__.lower().startswith("jerusalem")
