#! /usr/bin/env python
import argparse
from hera import datalayer
import json

parser = argparse.ArgumentParser()
parser.add_argument('command', nargs=1, type=str)
parser.add_argument('args', nargs='*', type=str)

args = parser.parse_args()


def load_handler(arguments):
    with open(arguments[4], 'r') as myFile:
        params = json.load(myFile)

    doc = dict(projectName=arguments[0],
               resource=arguments[1],
               dataFormat='string',
               type='LSM_template',
               desc=dict(params=params,
                         modelFolder=arguments[2],
                         version=arguments[3]
                         )
               )
    datalayer.Simulations.addDocument(**doc)


globals()['%s_handler' % args.command[0]](args.args)
