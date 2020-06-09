    def topography(self, infile, lowerleft1, lowerleft2, upperright1, upperright2, flat=None):
	dxdy = 150 # 10
	utls = geoPandas2STL.Topography(dxdy=float(dxdy))
	print("Converting...")
	stlstr, data = utls.Convert_shp_to_stl(shpfile=infile, solidname="Topography", flat=flat)
	return stlstr, np.min(data['Height']), np.max(data['Height']) 
	

    def buildings(self, outputfile, lowerleft1, lowerleft2, upperright1, upperright2, flat=None):

        maxheight = -500
	doc = FreeCAD.newDocument("Unnamed")
	shp = geopandas.read_file(outputfile)
	k=-1
	for j in range(len(shp)):  # converting al the buildings
	    try:
	    walls = shp['geometry'][j].exterior.xy
            except:
		continue
		print('bad j !!!', j)
	    if j % 100 == 0:
		print ('100j',j,k,len(shp)) # just to see that the function is still working

	    k = k + 1
	    wallsheight = shp['BLDG_HT'][j]
	    if flat is None:
	        altitude = shp['HT_LAND'][j] 
	    else:
	        altitude = flat
	    doc.addObject('Sketcher::SketchObject', 'Sketch' + str(j))
	    doc.Objects[2*k].Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, 0.000000), # 2*k-1
	                                              FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))

	    for i in range(len(walls[0])-1):
                doc.Objects[2*k].addGeometry(Part.Line(FreeCAD.Vector(walls[0][i], walls[1][i], altitude),
		                                      FreeCAD.Vector(walls[0][i+1], walls[1][i+1], altitude)))

	    doc.addObject("PartDesign::Pad", "Pad"+str(j))
	    doc.Objects[2 * k + 1].Sketch = doc.Objects[2*k]
	    buildingTopAltitude = wallsheight + altitude # wallsheight + altitude
	    maxheight = max(maxheight, buildingTopAltitude)
	    doc.Objects[2 * k + 1].Length = buildingTopAltitude # 35.000000
	    doc.Objects[2 * k + 1].Reversed = 0
	    doc.Objects[2 * k + 1].Midplane = 0
	    doc.Objects[2 * k + 1].Length2 = wallsheight # 100.000000
	    doc.Objects[2 * k + 1].Type = 0
	    doc.Objects[2 * k + 1].UpToFace = None
	    doc.recompute()  # maybe it can go outside the for loop

	doc.recompute()
	Mesh.export(doc.Objects, "file-buildings.stl")
 
        return maxheight


