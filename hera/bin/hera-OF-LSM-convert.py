#! /usr/bin/env python
import os
from hera.simulations.openFoam import OFObjectHome
import argparse

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('caseDir', nargs=1, type=str)
    parser.add_argument('outputdir', nargs=1, type=str)
    parser.add_argument('--workers', type=int,default=36,dest="workers")
    parser.add_argument('--overwrite', default=False, dest="overwrite",action="store_true")

    args = parser.parse_args()

    from dask.distributed import Client
    client = Client(n_workers=args.workers)

    data = OFObjectHome.loadLagrangianDataParallel(withVelocity=False,
                                                   withReleaseTimes=False,
                                                   withMass=True,
                                                   casePath=args.caseDir[0]).set_index("time").repartition(npartitions=20)


    outpath = os.path.abspath(args.outputdir[0])
    outName = os.path.join(outpath,"particles.parquet")

    if args.overwrite or not os.path.exists(outName):
        os.makedirs(outpath,exist_ok=True)
        data.to_parquet(outName)
    else:
        print("File exists, use overwrite to write new version")


    client.close()


