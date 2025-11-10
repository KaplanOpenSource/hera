# Tests for ToolkitHome.registerToolkit end-to-end:
# - registers a class on disk as a ToolkitDataSource
# - validates returned document fields
# - instantiates the class through DataHandler_Class (using the saved resource/desc)
# - verifies overwrite behavior

import os
import sys
import textwrap
from pathlib import Path
from importlib import import_module

from hera.toolkit import ToolkitHome, abstractToolkit, TOOLKIT_DATASOURCE_TYPE
from hera.datalayer.datahandler import DataHandler_Class


class _ToolkitHomeForTests(ToolkitHome, abstractToolkit):
    """
    Adapter that mixes ToolkitHome with abstractToolkit so we can reuse
    addDataSource/getMeasurementsDocuments without changing production code.
    """
    def __init__(self, project_name="UnitTestProject"):
        # abstractToolkit needs a toolkitName + projectName
        abstractToolkit.__init__(self, toolkitName="DynamicToolkits", projectName=project_name)
        # ToolkitHome seeds static toolkits
        ToolkitHome.__init__(self)


def _write_demo_toolkit(tmp_path: Path) -> tuple[str, str, str]:
    """
    Create a tiny toolkit on disk:
      <tmp_path>/mypkg/mytoolkit.py  with class DemoToolkit
    Returns:
      (resource_dir_to_add_to_sys_path, package_dir, classpath)
    """
    pkg = tmp_path / "mypkg"
    pkg.mkdir(parents=True, exist_ok=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "mytoolkit.py").write_text(textwrap.dedent("""
    class DemoToolkit:
        def __init__(self, projectName=None, filesDirectory=None, alpha=1):
            self.projectName = projectName
            self.filesDirectory = filesDirectory
            self.alpha = alpha

        def ping(self):
            return f"OK:{self.projectName}:{self.alpha}"
    """))
    resource_dir = str(tmp_path)                    # folder to add to sys.path for import
    package_dir = str(pkg)                          # sometimes stored as resource by registerToolkit
    classpath = "mypkg.mytoolkit.DemoToolkit"
    return resource_dir, package_dir, classpath


def _import_class(resource_dir: str, classpath: str):
    """
    Adds resource_dir to sys.path temporarily and imports the class from classpath.
    """
    if resource_dir not in sys.path:
        sys.path.insert(0, resource_dir)
    module_name, _, class_name = classpath.rpartition(".")
    mod = import_module(module_name)
    return getattr(mod, class_name)


def test_register_toolkit_creates_datasource_and_instantiates(tmp_path):
    th = _ToolkitHomeForTests(project_name="UnitTestProject")

    resource_dir, package_dir, classpath = _write_demo_toolkit(tmp_path)
    DemoToolkit = _import_class(resource_dir, classpath)

    # Register the class as a datasource (pass the class and projectName)
    doc = th.registerToolkit(
        toolkitclass=DemoToolkit,
        projectName=th.projectName,
        datasource_name="DemoToolkit_DS",
        params={"alpha": 7},
        version=(0, 0, 1),
        overwrite=True,  # tolerant across repeated runs
    )

    # Basic shape of the returned document
    assert hasattr(doc, "type") and doc.type == TOOLKIT_DATASOURCE_TYPE
    assert hasattr(doc, "resource") and os.path.isdir(doc.resource)
    assert hasattr(doc, "desc") and isinstance(doc.desc, dict)

    # Must include datasourceName, version, and classpath under desc
    assert doc.desc.get("datasourceName") == "DemoToolkit_DS"
    # version may be saved as list or tuple depending on backend
    v = doc.desc.get("version")
    assert tuple(v) == (0, 0, 1)
    assert "classpath" in doc.desc and doc.desc["classpath"].endswith(".DemoToolkit")

    # Instantiate the toolkit through DataHandler_Class using saved record
    obj = DataHandler_Class.getData(resource=doc.resource, desc=doc.desc)
    assert obj.__class__.__name__ == "DemoToolkit"
    assert obj.alpha == 7
    assert obj.ping().startswith("OK:")


def test_register_toolkit_overwrite_behavior(tmp_path):
    th = _ToolkitHomeForTests(project_name="UnitTestProject")
    resource_dir, package_dir, classpath = _write_demo_toolkit(tmp_path)
    DemoToolkit = _import_class(resource_dir, classpath)

    # First registration succeeds
    doc1 = th.registerToolkit(
        toolkitclass=DemoToolkit,
        projectName=th.projectName,
        datasource_name="DemoToolkit_DS",
        params={"alpha": 3},
        version=(0, 0, 1),
        overwrite=True,  # make the test idempotent
    )
    assert doc1 is not None

    # Second registration with the same name+version and overwrite=False should fail
    failed = False
    try:
        th.registerToolkit(
            toolkitclass=DemoToolkit,
            projectName=th.projectName,
            datasource_name="DemoToolkit_DS",
            params={"alpha": 4},
            version=(0, 0, 1),
            overwrite=False,
        )
    except Exception:
        failed = True
    assert failed, "Expected a failure on duplicate registration without overwrite=True"

    # But with overwrite=True it should succeed
    doc2 = th.registerToolkit(
        toolkitclass=DemoToolkit,
        projectName=th.projectName,
        datasource_name="DemoToolkit_DS",
        params={"alpha": 5},
        version=(0, 0, 1),
        overwrite=True,
    )
    assert doc2 is not None
