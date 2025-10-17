# hera/utils/data/cli_toolkit_repository.py
import argparse
import json
import os
import sys
from importlib import import_module

from hera.toolkit import ToolkitHome
from hera.utils.data.toolkit_repository import ToolkitRepository


# -----------------------------
# Helpers
# -----------------------------
def _str_to_version_tuple(s: str):
    """
    "0.0.1" -> (0, 0, 1)
    """
    parts = s.strip().split(".")
    try:
        return tuple(int(p) for p in parts)
    except Exception:
        raise ValueError(f"Invalid version string '{s}'. Expected format like 0.0.1")


def _ensure_on_syspath(path: str):
    """Add directory to sys.path if needed."""
    if not path:
        return
    abspath = os.path.abspath(path)
    if os.path.isdir(abspath) and abspath not in sys.path:
        sys.path.append(abspath)


def _import_class(classpath: str):
    """
    Import a class given a fully-qualified class path, e.g. 'pkg.module.ClassName'.
    """
    mod_name, _, cls_name = classpath.rpartition(".")
    if not mod_name or not cls_name:
        raise ValueError(
            f"Invalid --classpath '{classpath}'. Expected something like 'pkg.module.ClassName'."
        )
    mod = import_module(mod_name)
    return getattr(mod, cls_name)


def _datasource_exists(project: str, repository: str, name: str) -> bool:
    """
    Check if a ToolkitDataSource document exists for (project, repository, name).
    """
    try:
        from hera.datalayer.project import Project
        proj = Project(projectName=project)
        docs = proj.getMeasurementsDocuments(
            type="ToolkitDataSource",
            repository=repository,
            datasourceName=name,
        )
        return bool(docs)
    except Exception:
        # If the project layer is unavailable for some reason, be safe and say False
        return False


# -----------------------------
# Commands
# -----------------------------
def cmd_register_datasource(args: argparse.Namespace) -> None:
    """
    Register a toolkit as a ToolkitDataSource by inspecting a class object.
    Uses ToolkitHome.registerToolkit() and stores: resource (folder), classpath, parameters, repository, version.
    """
    # Optional: put resource on sys.path so import can succeed
    if args.resource:
        _ensure_on_syspath(args.resource)

    # Parse params JSON (optional)
    params = json.loads(args.params) if args.params else {}

    # If user asked to ignore when exists, check first
    if args.ignore and _datasource_exists(args.project, args.repository, args.name):
        print(f"(ignored) Datasource '{args.name}' already exists in repository '{args.repository}' for project '{args.project}'.")
        return

    # Import the class object from --classpath (e.g. mypkg.mymod.MyClass)
    toolkit_cls = _import_class(args.classpath)

    # Register
    th = ToolkitHome()
    doc = th.registerToolkit(
        toolkitclass=toolkit_cls,
        datasource_name=args.name,
        params=params,
        version=_str_to_version_tuple(args.version),
        overwrite=args.overwrite,
        projectName=args.project,
        repositoryName=args.repository,  # <<< important: persist repository in the document
    )

    # Print summary
    print("Registered datasource:")
    print("  project     :", args.project)
    print("  repository  :", args.repository)
    print("  name        :", args.name)
    print("  classpath   :", doc.desc.get("classpath"))
    print("  resource    :", doc.resource)
    print("  parameters  :", doc.desc.get("parameters"))
    print("  version     :", doc.desc.get("version"))


def cmd_add_doc(args: argparse.Namespace) -> None:
    """
    Add/overwrite a 'Toolkit' document in the DB (registry entry),
    using the ToolkitRepository.
    """
    repo = ToolkitRepository(args.project)
    _ = repo.addToolkitDocument(
        toolkitName=args.name,
        cls=args.cls,
        source=args.source,
        tk_type=args.type,
        description=args.description,
        overwrite=args.overwrite,
    )
    print("Added Toolkit document:")
    print("  project     :", args.project)
    print("  toolkit     :", args.name)
    print("  cls         :", args.cls)
    print("  source      :", args.source)
    print("  type        :", args.type)
    print("  description :", args.description)


def cmd_print(args: argparse.Namespace) -> None:
    """
    Print the toolkits table for a project (static + dynamic).
    """
    th = ToolkitHome()
    df = th.getToolkitTable(args.project)
    if args.filter_source:
        df = df[df["source"] == args.filter_source]
    if args.filter_type:
        df = df[df["type"] == args.filter_type]

    if args.format == "csv":
        print(df.to_csv(index=False))
    else:
        if df.empty:
            print("(no toolkits found)")
        else:
            try:
                from tabulate import tabulate  # optional
                print(tabulate(df, headers="keys", tablefmt="github", showindex=False))
            except Exception:
                print(df.to_string(index=False))


def cmd_set_default_repository(args):
    th = ToolkitHome()
    th.setDefaultRepository(projectName=args.project, repositoryName=args.repository, overwrite=True)
    print(f"✔ defaultRepository for project '{args.project}' set to '{args.repository}'")

def cmd_get_datasource(args):
    th = ToolkitHome()
    version = None
    if args.version:
        version = tuple(int(x) for x in args.version.split("."))
    doc = th.getDatasourceDocument(
        projectName=args.project,
        datasourceName=args.name,
        repositoryName=args.repository,  # may be None -> fallback to default
        version=version,
    )
    if not doc:
        print("(!) Datasource not found.")
        return
    print("Found datasource:")
    print("  project     :", args.project)
    print("  repository  :", doc.desc.get("repository"))
    print("  name        :", doc.desc.get("datasourceName"))
    print("  classpath   :", doc.desc.get("classpath"))
    print("  version     :", doc.desc.get("version"))
    print("  resource    :", doc.resource)
    print("  parameters  :", doc.desc.get("parameters"))


# -----------------------------
# CLI
# -----------------------------
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
    r.add_argument("--overwrite", action="store_true", help="Overwrite if datasource already exists")
    r.add_argument("--ignore", action="store_true", help="If exists, do nothing (no-op)")
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

    # set-default-repository
    sdr = sub.add_parser("set-default-repository", help="Set the project's default repository")
    sdr.add_argument("--project", "-p", required=True, help="Project name")
    sdr.add_argument("--repository", "-r", required=True, help="Default repository name to set")
    sdr.set_defaults(func=cmd_set_default_repository)

    # get-datasource
    gd = sub.add_parser("get-datasource", help="Fetch a ToolkitDataSource by (repository?, name[, version])")
    gd.add_argument("--project", "-p", required=True, help="Project name")
    gd.add_argument("--name", "-n", required=True, help="Datasource name")
    gd.add_argument("--repository", "-r", required=False, help="Repository name (optional; falls back to default)")
    gd.add_argument("--version", "-v", required=False, help="Version like 0.0.1 (optional)")
    gd.set_defaults(func=cmd_get_datasource)

    # print
    ls = sub.add_parser("print", help="Print toolkits table (static + dynamic)")
    ls.add_argument("--project", "-p", required=True, help="Project name")
    ls.add_argument("--filter-source", choices=["internal", "experiment", "external"], help="Filter by source")
    ls.add_argument("--filter-type", choices=["measurements", "simulations"], help="Filter by type")
    ls.add_argument("--format", choices=["table", "csv"], default="table", help="Output format")
    ls.set_defaults(func=cmd_print)

    return p


def main(argv=None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
