#! /usr/bin/env python
import argparse
from hera import datalayer

parser = argparse.ArgumentParser()
parser.add_argument('command', nargs=1, type=str)
parser.add_argument('args', nargs='*', type=str)

args = parser.parse_args()


def load_handler(arguments):
    extents = eval("%s" % arguments[1])
    doc = dict(projectName=arguments[0],
               resource=arguments[3],
               dataFormat='image',
               type='GIS',
               desc=dict(locationName=arguments[2],
                         left=extents[0],
                         right=extents[1],
                         bottom=extents[2],
                         top=extents[3]
                         )
               )
    datalayer.Measurements.addDocument(**doc)


globals()['%s_handler' % args.command[0]](args.args)
