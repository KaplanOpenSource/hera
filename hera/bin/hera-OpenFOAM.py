import argparse
import json

import pandas
import os
import shutil
import sys
version = sys.version_info[0]
if version == 2:
    from hera import openfoam
else:
    from hera import datalayer
    import dask.dataframe as dd

parser = argparse.ArgumentParser()
parser.add_argument('command', nargs=1, type=str)
parser.add_argument('args', nargs='*', type=str)
parser.add_argument('-keepHDF', action='store_false')

args = parser.parse_args()


commands = ["load"]

def load_handler(arguments):
    """
    Converts hdf files with OpenFOAM results to parquet format, adds the paths to the database.
    Unless -keepHDF is used, deletes the hdf files.
    path is the path of the directory in which the hdf directory rests.
    name is the name used for the formation of the hdf files.
    projectName is used for saving the data to the database.
    """

    path = arguments[0]
    name = arguments[1]
    projectName = arguments[2]

    keys = pandas.HDFStore("%s/%s/hdf/%s_0.hdf" % (path, name, name)).keys()
    fixed_keys = []
    for key in keys:
        if key[0]=="/":
            fixed_keys.append(key[1:])
        else:
            fixed_keys.append(key)

    if not os.path.isdir("%s/%s/parquet" % (path, name)):
        os.makedirs("%s/%s/parquet" % (path, name))

    print(fixed_keys)
    with open("%s/%s/meta.json" % (path, name)) as json_file:
        jsondata = json.load(json_file)
    metadata = jsondata["metadata"] # Reading the metadata

    # Making a list of the pipeline filter names

    names = {}

    def makenames(pipeline, former):
        for pipe in pipeline:
            names.update({pipe: former+"_"+pipe})
            print(pipe)
            if "downstream" in pipeline[pipe].keys():
                makenames(pipeline[pipe]["downstream"], former=former+pipe)
    makenames(jsondata["pipelines"],"")


    # Making the parquet files

    for key in fixed_keys:
        data = dd.read_hdf("%s/%s/hdf/%s_0.hdf" % (path, name, name), key=key)
        for i in range(1, len([name for name in os.listdir('%s/%s/hdf' % (path, name)) if os.path.isfile(os.path.join('%s/%s/hdf' % (path, name), name))])):
            new_data = dd.read_hdf("%s/%s/hdf/%s_%s.hdf" % (path, name, name, (i)), key=key)
            data = dd.concat([data, new_data],interleave_partitions=True).reset_index().set_index("time").repartition(partition_size="100Mb")
        data.to_parquet("%s/%s/parquet/%s.parquet" % (path, name, key))

        # Adding a link to the parquet to the database

        datalayer.Measurements.addDocument(projectName=projectName, desc=dict(filter=key, pipeline=names[key], **metadata), type="OFsimulation",
                                           resource="%s/%s/parquet/%s.parquet" % (path, name, key), dataFormat="parquet")

    #   Delete hdf directory
    if args.keepHDF:
        shutil.rmtree("%s/%s/hdf" % (path, name), ignore_errors=True)
        shutil.rmtree("%s/%s/meta.json" % (path, name), ignore_errors=True)

def executePipeline_handler(arguments):

    JSONpath = arguments[0]
    name = arguments[1]
    casePath = arguments[2]
    caseType = arguments[3] if len(arguments) > 3 else "Decomposed Case"
    servername = arguments[4] if len(arguments) > 4 else None

    with open(JSONpath) as json_file:
        data = json.load(json_file)
    vtkpipe = openfoam.VTKpipeline(name=name, pipelineJSON=data, casePath=casePath, caseType=caseType, servername=servername)
    vtkpipe.execute("mainReader")

    # (name
    #  ="test", pipelineJSON=data, casePath="/home/ofir/Projects/openFoamUsage/askervein", caseType="Reconstructed Case")

globals()['%s_handler' % args.command[0]](args.args)