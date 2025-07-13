#!/home/ilay/hera/heraenv/bin/python3.9

from PyFoam.Applications.ChangePython import changePython

changePython("pvpython","PVSnapshot",options=["--native"])
