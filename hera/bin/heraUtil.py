import argparse
import json
from hera import datalayer

Measurements = datalayer.Measurements
Simulations = datalayer.Simulations
Analysis = datalayer.Cache

parser = argparse.ArgumentParser()
parser.add_argument('--dump', dest='dump', action='store_true', help='')
parser.add_argument('--outputFile', dest='outputFile', default=None, help='')
parser.add_argument('--projectName', dest='projectName', default=None, help='')
parser.add_argument('--type', dest='type', default=None, help='')
parser.add_argument('--queryFile', dest='queryFile', default=None, help='')
parser.add_argument('--load', dest='load', action="store_true", help='')
parser.add_argument('--dataFile', dest='dataFile', default=None, help='')

args = parser.parse_args()

if args.dump:
    if args.outputFile is None:
        parser.error("Must use --outputFile with --dump")
    else:
        if args.queryFile is None:
            if args.projectName is None:
                parser.error("Must use --projectName with --dump")
            elif args.type is None:
                parser.error("Must use --type when using --projectName with --dump")
            elif args.type not in ['Measurements', 'Simulations', 'Analysis']:
                parser.error("--type must get one of the following values: 'Measurements', 'Simulations', 'Analysis'")
            else:
                with open(args.outputFile, 'w') as myFile:
                    json.dump(globals()[args.type].getDocumentsAsDict(args.projectName), myFile, indent=4, sort_keys=True)
        else:
            with open(args.queryFile, 'r') as myFile:
                query = json.load(myFile)
            if 'projectName' not in query.keys():
                raise ValueError("'projectName' key must be given in the query file")
            elif 'type' not in query.keys():
                raise ValueError("'type' key must be given in the query file")
            else:
                with open(args.outputFile, 'w') as myFile:
                    json.dump(globals()[query.pop('type')].getDocumentsAsDict(**query), myFile, indent=4, sort_keys=True)

elif args.load:
    if args.dataFile is None:
        parser.error("Must use --dataFile with --load")
    else:
        with open(args.dataFile, 'r') as myFile:
            docList = json.load(myFile)['documents']
        for doc in docList:
            docCls = doc.pop('_cls').split('.')[-1]
            globals()[docCls].addDocumentFromJSON(json.dumps(doc))