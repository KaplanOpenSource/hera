import os
import sys
import pytest


def _import_datahandler():
    """
    Small helper:
    - Adds the root of the Hera repository to sys.path
    - Imports and returns DataHandler_Class from the data layer
    - Also returns the repo root path for convenience.
    """
    # __file__ points to: hera/tests/dynamic_loading_tests_pack/test_*.py
    # going three levels up => /path/to/hera
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))

    # Ensure Hera root is in sys.path for dynamic imports
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    # Import the dynamic loading handler from Hera
    from hera.datalayer.datahandler import DataHandler_Class
    return DataHandler_Class, repo_root


@pytest.mark.integration
def test_dynamic_toolkit_datasource_can_be_loaded():
    """
    Full integration test for dynamic toolkit loading:

    Goal:
    Validate that Hera can dynamically load a toolkit class from source code,
    using only a filesystem resource path + classpath string.

    Steps:
    1. Ensure the dummy toolkit package exists under:
       hera/tests/dynamic_loading/dummy_toolkit
    2. Call DataHandler_Class.getData with:
       - resource = path to the package folder
       - desc["classpath"] = fully qualified class path
    3. Verify that:
       - The class is dynamically imported
       - An instance is created
       - It is of type DummyToolkit
       - The ping() method works
    """

    DataHandler_Class, repo_root = _import_datahandler()

    # Full path to the dummy toolkit package inside the Hera repository
    dummy_pkg_path = os.path.join(
        repo_root,
        "hera",
        "tests",
        "dynamic_loading",
        "dummy_toolkit",
    )

    # Ensure the dummy toolkit directory exists
    assert os.path.isdir(dummy_pkg_path), (
        f"Dummy toolkit package not found at {dummy_pkg_path}"
    )

    # Fully qualified classpath to import dynamically
    classpath = "hera.tests.dynamic_loading.dummy_toolkit.DummyToolkit"

    # Perform dynamic loading using DataHandler_Class.getData
    obj = DataHandler_Class.getData(
        resource=dummy_pkg_path,     # points to the package folder
        desc={"classpath": classpath},  # tells Hera which class to load
    )

    # Assertions to validate dynamic loading success
    assert obj is not None, "Expected a toolkit instance, got None."
    assert obj.__class__.__name__ == "DummyToolkit"
    assert hasattr(obj, "ping"), "DummyToolkit instance has no ping() method."
    assert obj.ping().startswith("dummy-toolkit-ok:"), (
        "ping() did not return the expected output."
    )
