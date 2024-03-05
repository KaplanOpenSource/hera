#! /usr/bin/env python


if __name__ == "__main__":
    conf = loadJSON(os.path.join(pathlib.Path(__file__).parent, "datasources.json"))

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subparser_name")

    load_parser = subparsers.add_parser('load')
    load_parser.add_argument('projectName', type=str, help='The project to load the data')
    load_parser.add_argument('dataType', type=str, help='The data type to load')
    load_parser.add_argument('source', type=str, help='The name of the datasource to load')
    load_parser.set_defaults(func=load)

    list_parser = subparsers.add_parser('list')
    list_parser.add_argument('dataType', type=str, default=None, nargs='*')
    load_parser.set_defaults(func=list)

    loadAll_parser = subparsers.add_parser('loadAll')
    loadAll_parser.add_argument('projectName', default=None, type=str, help='The project to load the data')
    loadAll_parser.add_argument('--overwrite', dest="overwrite", action="store_true", default=False,
                                help='Forces overwriting of the datasource, it exists')
    loadAll_parser.set_defaults(func=loadAll)

    loadAll_parser = subparsers.add_parser('loadDefaultProject')
    loadAll_parser.set_defaults(overwrite=True)
    loadAll_parser.set_defaults(projectName=None)
    loadAll_parser.set_defaults(func=loadAll)

    args = parser.parse_args()
    args.func(args, conf)


