# hera/utils/data/cli_toolkit_repository.py
# -----------------------------------------------------------------------------
# Generic CLI for dynamic Toolkit registration and inspection.
# This version:
#   1) Adds --repository (-R) to `register-datasource` and passes it into
#      ToolkitHome.registerToolkit(..., repositoryName=...).
#   2) Supports project-level default repository management:
#        - default-repo set/get
#      If --repository is omitted at registration time, the CLI falls back
#      to the project's default repository (and errors if none is set).
#   3) Robust class import from --classpath, optional --resource to extend
#      sys.path on the fly, and JSON parsing for --params.
#   4) Keeps the existing `add-doc` (static registry entry) and `print`
#      (static + dynamic toolkits table).
# -----------------------------------------------------------------------------

import argparse
import json
import os
import sys
from importlib import import_module
from typing import Any, Dict, Tuple

from hera.toolkit import ToolkitHome
from hera.utils.data.toolkit_repository import ToolkitRepository


# ------------------------------- Utilities -----------------------------------

def _str_to_version_tuple(s: str) -> Tuple[int, ...]:
    """Convert a dotted string (e.g., '0.0.1') into a version tuple, e.g., (0, 0, 1)."""
    parts = s.strip().split(".")
    try:
        return tuple(int(p) for p in parts if p != "")
    except Exception as exc:
        raise ValueError(
            f"Invalid version string '{s}'. Expected format like '0.0.1'"
        ) from exc


def _parse_params_json(params_str: str) -> Dict[str, Any]:
    """Parse JSON string for constructor params. Returns {} for empty/whitespace."""
    params_str = (params_str or "").strip()
    if not params_str:
        return {}
    try:
        obj = json.loads(params_str)
    except Exception as exc:
        raise ValueError(
            f"--params must be a valid JSON object, got: {params_str}"
        ) from exc
    if not isinstance(obj, dict):
        raise ValueError("--params must be a JSON object (e.g. '{\"alpha\": 7}')")
    return obj


def _import_toolkit_class(classpath: str, resource: str = "") -> type:
    """
    Import a class by fully-qualified classpath like 'pkg.mod.ClassName'.
    Optionally extends sys.path with --resource before importing.
    """
    classpath = (classpath or "").strip()
    if not classpath or "." not in classpath:
        raise ValueError(
            f"--classpath must be a fully-qualified path like 'pkg.mod.ClassName', got: {classpath!r}"
        )

    if resource:
        # Allow ad-hoc injection of a local folder that contains the module.
        resource = os.path.abspath(resource)
        if os.path.isdir(resource) and resource not in sys.path:
            sys.path.insert(0, resource)

    mod_name, _, cls_name = classpath.rpartition(".")
    if not mod_name or not cls_name:
        raise ValueError(
            f"--classpath must end with a class name (e.g. 'pkg.mod.ClassName'), got: {classpath!r}"
        )

    try:
        mod = import_module(mod_name)
    except Exception as exc:
        raise ImportError(
            f"Failed to import module '{mod_name}' for classpath '{classpath}'."
        ) from exc

    try:
        cls = getattr(mod, cls_name)
    except AttributeError as exc:
        raise ImportError(
            f"Module '{mod_name}' does not define class '{cls_name}' for classpath '{classpath}'."
        ) from exc

    if not isinstance(cls, type):
        raise TypeError(f"Resolved object '{cls_name}' is not a class (classpath={classpath}).")

    return cls


# ------------------------------- Commands ------------------------------------

def cmd_register_datasource(args: argparse.Namespace) -> None:
    """
    Register a dynamic ToolkitDataSource into the project's repository via ToolkitHome.
    This persists a DB-backed 'ToolkitDataSource' record (dynamic registry).
    """
    # Import toolkit class (and extend sys.path if --resource was provided)
    toolkit_cls = _import_toolkit_class(args.classpath, args.resource)

    # Parse constructor parameters
    params = _parse_params_json(args.params)

    # If user asked to ignore when exists, check first
    if args.ignore and _datasource_exists(args.project, args.repository, args.name):
        print(f"(ignored) Datasource '{args.name}' already exists in repository '{args.repository}' for project '{args.project}'.")
        return

    # Import the class object from --classpath (e.g. mypkg.mymod.MyClass)
    toolkit_cls = _import_class(args.classpath)

    # Register
    th = ToolkitHome()

    # Resolve repository: CLI flag takes precedence; otherwise use project default
    repo = (args.repository or "").strip()
    if not repo:
        repo = th.getDefaultRepository(projectName=args.project)
        if not repo:
            raise ValueError(
                "No --repository provided and no default repository set for this project. "
                "Set a default repository first:\n"
                "  python -m hera.utils.data.cli_toolkit_repository default-repo set "
                f"-p {args.project} -R <REPOSITORY_NAME>"
            )

    # Perform registration (this requires repositoryName)
    doc = th.registerToolkit(
        toolkitclass=toolkit_cls,
        datasource_name=args.name,
        params=params,
        version=_str_to_version_tuple(args.version),
        overwrite=bool(args.overwrite),
        projectName=args.project,
        repositoryName=repo,
    )

    # Human-readable summary
    print("Registered datasource:")
    print("  project     :", args.project)
    print("  repository  :", repo)
    print("  name        :", args.name)
    print("  classpath   :", doc.desc.get("classpath"))
    print("  resource    :", doc.resource)
    print("  parameters  :", doc.desc.get("parameters"))
    print("  version     :", doc.desc.get("version"))


def cmd_add_doc(args: argparse.Namespace) -> None:
    """
    Create/overwrite a 'ToolkitDataSource' measurements document so that the
    toolkit appears in the unified toolkits table, without requiring an importable class.
    If --cls is provided but cannot be imported, we silently omit 'classpath' to keep
    ToolkitHome.getToolkit() safe (it will fall back to the shim).
    """
    # Resolve safe classpath (optional)
    classpath = (args.cls or "").strip()
    keep_classpath = False
    if classpath:
        mod_name, _, cls_name = classpath.rpartition(".")
        if mod_name and cls_name:
            try:
                # Try importing to verify the classpath really exists
                mod = import_module(mod_name)
                getattr(mod, cls_name)
                keep_classpath = True
            except Exception:
                # Do not raise; we will omit classpath so getToolkit falls back to the shim
                keep_classpath = False

    # Build the measurements document descriptor
    desc = {
        "datasourceName": args.name,
        "toolkit": args.name,
        "toolkitType": args.type,       # keeps alignment with getToolkitTable()
        "source": args.source,          # "internal" / "experiment" / "external"
        "description": args.description or "",
    }
    if keep_classpath:
        desc["classpath"] = classpath   # only if verified importable

    # Upsert via Project API
    repo = ToolkitRepository(args.project)               # helper for lookups
    existing = repo.getToolkitDocument(args.name)        # returns measurements doc or None

    if existing and not args.overwrite:
        print(f"(exists) Toolkit '{args.name}' already present; use --overwrite to replace")
        return

    # If exists and overwrite requested -> delete old doc
    if existing:
        try:
            existing.delete()
        except Exception:
            pass

    # Insert a fresh measurements document with type="ToolkitDataSource"
    proj = repo._project  # underlying Project
    proj.addMeasurementsDocument(
        type="ToolkitDataSource",
        dataFormat="Class",
        desc=desc,
    )

    print("Added Toolkit document:")
    print("  project     :", args.project)
    print("  toolkit     :", args.name)
    print("  cls(saved)  :", desc.get("classpath", "(omitted)"))
    print("  source      :", args.source)
    print("  type        :", args.type)
    print("  description :", args.description or "")


def cmd_print(args: argparse.Namespace) -> None:
    """
    Print the toolkits table for a project (static + dynamic union).
    Optional filters by 'source' and 'type'. Supports CSV or a simple table.
    """
    th = ToolkitHome()
    df = th.getToolkitTable(args.project)

    if args.filter_source:
        df = df[df["source"] == args.filter_source]
    if args.filter_type:
        df = df[df["type"] == args.filter_type]

    if args.format == "csv":
        print(df.to_csv(index=False))
        return

    if df.empty:
        print("(no toolkits found)")
        return

    # Pretty textual table using tabulate (if available); otherwise fallback to DataFrame printing.
    try:
        from tabulate import tabulate  # optional dependency
        print(tabulate(df, headers="keys", tablefmt="github", showindex=False))
    except Exception:
        print(df.to_string(index=False))


def cmd_default_repo_set(args: argparse.Namespace) -> None:
    """
    Set the project's default repository. Future registrations can omit --repository.
    """
    th = ToolkitHome()
    th.setDefaultRepository(projectName=args.project, repositoryName=args.repository, overwrite=True)
    print(f"Default repository for project '{args.project}' set to '{args.repository}'")


def cmd_default_repo_get(args: argparse.Namespace) -> None:
    """
    Get the project's default repository (if any).
    """
    th = ToolkitHome()
    repo = th.getDefaultRepository(projectName=args.project)
    print(repo if repo else "(no default repository set)")


# --------------------------------- Parser ------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hera-repository (python -m hera.utils.data.cli_toolkit_repository)",
        description="Register/inspect Hera toolkits dynamically."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # register-datasource
    r = sub.add_parser("register-datasource", help="Register a toolkit class as a ToolkitDataSource")
    r.add_argument("--project", "-p", required=True, help="Project name")
    r.add_argument("--repository", "-R", required=True, help="Repository name to store the datasource")  # <<< required
    r.add_argument("--name", "-n", required=True, help="Datasource name to register under (unique within repository)")
    r.add_argument("--classpath", "-c", required=True, help="Fully-qualified class path (e.g. mypkg.mymod.MyClass)")
    r.add_argument("--resource", "-r", default="", help="Optional folder to add to sys.path before import")
    r.add_argument("--params", "-P", default="", help='JSON dict of constructor kwargs (e.g. \'{"alpha": 7}\')')
    r.add_argument("--version", "-v", default="0.0.1", help="Version, e.g. 0.0.1")
    r.add_argument("--overwrite", action="store_true", help="Overwrite if exists")
    r.add_argument("--repository", "-R", default="", help="Repository name; if omitted, uses project's default repository")
    r.set_defaults(func=cmd_register_datasource)

    # add-doc (static registry entry via ToolkitRepository)
    a = sub.add_parser("add-doc", help="Add/overwrite a Toolkit registry document")
    a.add_argument("--project", "-p", required=True, help="Project name")
    a.add_argument("--name", "-n", required=True, help="Toolkit name (registry key)")
    a.add_argument("--cls", required=True, help="Python class path, e.g. hera.simulations.gaussian.toolkit.gaussianToolkit")
    a.add_argument("--source", required=True, choices=["internal", "experiment", "external"], help="Toolkit source")
    a.add_argument("--type", required=True, choices=["measurements", "simulations"], help="Toolkit type")
    a.add_argument("--description", default="", help="Optional description")
    a.add_argument("--overwrite", action="store_true", help="Overwrite if exists")
    a.set_defaults(func=cmd_add_doc)

    # print (unified static + dynamic view)
    ls = sub.add_parser("print", help="Print toolkits table (static + dynamic)")
    ls.add_argument("--project", "-p", required=True, help="Project name")
    ls.add_argument("--filter-source", choices=["internal", "experiment", "external"], help="Filter by source")
    ls.add_argument("--filter-type", choices=["measurements", "simulations"], help="Filter by type")
    ls.add_argument("--format", choices=["table", "csv"], default="table", help="Output format")
    ls.set_defaults(func=cmd_print)

    # default-repo (project-level default repository management)
    dr = sub.add_parser("default-repo", help="Manage default repository for a project")
    dr_sub = dr.add_subparsers(dest="op", required=True)

    dr_set = dr_sub.add_parser("set", help="Set default repository")
    dr_set.add_argument("--project", "-p", required=True, help="Project name")
    dr_set.add_argument("--repository", "-R", required=True, help="Repository name to set as default")
    dr_set.set_defaults(func=cmd_default_repo_set)

    dr_get = dr_sub.add_parser("get", help="Get default repository")
    dr_get.add_argument("--project", "-p", required=True, help="Project name")
    dr_get.set_defaults(func=cmd_default_repo_get)

    return p


# ---------------------------------- Main -------------------------------------

def main(argv=None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
