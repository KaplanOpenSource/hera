# hera/utils/data/cli_toolkit_repository.py
import argparse
import json
import os
import sys
from importlib import import_module

from hera.toolkit import ToolkitHome
from hera.utils.data.toolkit_repository import ToolkitRepository


def _str_to_version_tuple(s: str):
    # "0.0.1" -> (0, 0, 1)
    parts = s.strip().split(".")
    try:
        return tuple(int(p) for p in parts)
    except Exception:
        raise ValueError(f"Invalid version string '{s}'. Expected format like 0.0.1")


def cmd_register_datasource(args: argparse.Namespace) -> None:
    """
    Register a toolkit as a ToolkitDataSource by inspecting a class object.
    This uses ToolkitHome.registerToolkit, which:
      - figures out the module file path via inspect
      - stores resource (folder), classpath and parameters
    """
    # Optionally add a folder to sys.path so the import can work
    if args.resource:
        abs_path = os.path.abspath(args.resource)
        if os.path.isdir(abs_path) and abs_path not in sys.path:
            sys.path.append(abs_path)

    # Import the class object from --classpath (e.g. mypkg.mymod.MyClass)
    mod_name, _, cls_name = args.classpath.rpartition(".")
    if not mod_name or not cls_name:
        raise ValueError(
            f"Invalid --classpath '{args.classpath}'. Expected something like 'pkg.module.ClassName'."
        )
    mod = import_module(mod_name)
    toolkit_cls = getattr(mod, cls_name)

    params = json.loads(args.params) if args.params else {}

    th = ToolkitHome()
    doc = th.registerToolkit(
        toolkitclass=toolkit_cls,
        datasource_name=args.name,
        params=params,
        version=_str_to_version_tuple(args.version),
        overwrite=args.overwrite,
        projectName=args.project,
    )

    print("Registered datasource:")
    print("  project     :", args.project)
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
    doc = repo.addToolkitDocument(
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
        # Pretty-ish text output
        if df.empty:
            print("(no toolkits found)")
        else:
            from tabulate import tabulate  # optional dependency; if missing, fallback to print(df)
            try:
                print(tabulate(df, headers="keys", tablefmt="github", showindex=False))
            except Exception:
                print(df.to_string(index=False))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hera-repository (python -m hera.utils.data.cli_toolkit_repository)",
        description="Register/inspect Hera toolkits dynamically."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # register-datasource (dynamic ToolkitDataSource via ToolkitHome.registerToolkit)
    r = sub.add_parser("register-datasource", help="Register a toolkit class as a ToolkitDataSource")
    r.add_argument("--project", "-p", required=True, help="Project name")
    r.add_argument("--name", "-n", required=True, help="Datasource name to register under")
    r.add_argument("--classpath", "-c", required=True, help="Fully-qualified class path (e.g. mypkg.mymod.MyClass)")
    r.add_argument("--resource", "-r", default="", help="Optional folder to add to sys.path before import")
    r.add_argument("--params", "-P", default="", help='JSON dict of constructor kwargs (e.g. \'{"alpha": 7}\')')
    r.add_argument("--version", "-v", default="0.0.1", help="Version, e.g. 0.0.1")
    r.add_argument("--overwrite", action="store_true", help="Overwrite if exists")
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
