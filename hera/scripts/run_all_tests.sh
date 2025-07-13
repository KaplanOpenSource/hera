#!/bin/bash

echo "🚀 Running all Hera unit tests..."

# Low Frequency Meteorology
echo -e "\n🧪 Running test_unit_lowfreq_toolkit.py..."
python3 -m measurements.meteorology.lowfreqdata.test_unit_lowfreq_toolkit

# High Frequency Meteorology
echo -e "\n🧪 Running test_unit_highfreq_toolkit.py..."
python3 -m measurements.meteorology.highfreqdata.test_unit_highfreq_toolkit

# GIS Raster - Landcover
echo -e "\n🧪 Running test_unit_test_gis_raster_landcover.py..."
python3 -m measurements.GIS.raster.test_unit_test_gis_raster_landcover

# GIS Raster - Topography
echo -e "\n🧪 Running test_unit_test_gis_raster_topography.py..."
python3 -m measurements.GIS.raster.test_unit_test_gis_raster_topography

# GIS Vector - Demography
echo -e "\n🧪 Running test_unit_demography_toolkit.py..."
python3 -m measurements.GIS.vector.test_unit_demography_toolkit

echo -e "\n✅ All tests finished."
