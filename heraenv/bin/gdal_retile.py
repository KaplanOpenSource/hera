#!/home/ilay/hera/heraenv/bin/python3.9

import sys

from osgeo.gdal import UseExceptions, deprecation_warn

# import osgeo_utils.gdal_retile as a convenience to use as a script
from osgeo_utils.gdal_retile import *  # noqa
from osgeo_utils.gdal_retile import main

UseExceptions()

deprecation_warn("gdal_retile")
sys.exit(main(sys.argv))
