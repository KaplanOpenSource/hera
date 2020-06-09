import argparse
import json
from hera.datalayer import Simulations

parser = argparse.ArgumentParser()
parser.add_argument('--jsonFile', dest='jsonFile', required=True, help='The document JSON file to upload')

args = parser.parse_args()

with open(args.jsonFile, 'r') as myFile:
    doc = json.load(myFile)

Simulations.addDocument(**doc)