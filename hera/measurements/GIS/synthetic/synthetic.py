import FreeCAD
import Mesh

# Should run with "source initana2.sh"

class synthetic():

    def createSynthetic(self, stlpath, domainx, domainy, buildingx, buildingy, buildingz, gap):
        """
        Create a region with homogenous buildings, it is easier to check the code with synthetic environment
        
        param stlpath - the path where we save the stl
        param domainx - the width of the domain in the x direction
        param domainy - the width of the domain in the x direction
        param buildingx - the width of each building in the x direction
        param buildingy - the width of each building in the y direction
        param buildingz - the height of each building
        param gap - the distance between the buildings in the x direction and the y direction, the domain will start with a gap
        
        return:
        """
        doc = FreeCAD.newDocument("Unnamed")
        print('create synthetic')
        nx = 0
        km=-1
        while nx < domainx - gap - buildingx:
#            nx += gap 
            ny = 0                
            while ny < domainy - gap - buildingy:
#                ny += gap
                km += 1
                doc.addObject('Sketcher::SketchObject', 'Sketch' + str(km))
                doc.Objects[2*km].Placement = FreeCAD.Placement(FreeCAD.Vector(0.000000, 0.000000, 0.000000), # 2*k-1
		                                                           FreeCAD.Rotation(0.000000, 0.000000, 0.000000, 1.000000))
                doc.Objects[2*km].addGeometry(Part.Line(FreeCAD.Vector(nx, ny, 0),
		                                                  FreeCAD.Vector(nx, ny+buildingy, 0)))
                doc.Objects[2*km].addGeometry(Part.Line(FreeCAD.Vector(nx, ny+buildingy, 0),
		                                                  FreeCAD.Vector(nx+buildingx, ny+buildingy, 0)))
                doc.Objects[2*km].addGeometry(Part.Line(FreeCAD.Vector(nx+buildingx, ny+buildingy, 0),
		                                                  FreeCAD.Vector(nx+buildingx, ny, 0)))
                doc.Objects[2*km].addGeometry(Part.Line(FreeCAD.Vector(nx+buildingx, ny, 0),
		                                                  FreeCAD.Vector(nx, ny, 0)))
		                            
                doc.addObject("PartDesign::Pad", "Pad"+str(km))
                doc.Objects[2 * km + 1].Sketch = doc.Objects[2*km]
                doc.Objects[2 * km + 1].Length = buildingz # 35.000000
                doc.Objects[2 * km + 1].Reversed = 0
                doc.Objects[2 * km + 1].Midplane = 0
                doc.Objects[2 * km + 1].Length2 = buildingz # 100.000000
                doc.Objects[2 * km + 1].Type = 0
                doc.Objects[2 * km + 1].UpToFace = None
                doc.recompute()
                ny += buildingy + gap
            nx += buildingx + gap
		                
        doc.recompute()
		                # Mesh.export([doc.getObject("Extrude")], u"/home/nirb/regions/menuO1.stl")
        Mesh.export(doc.Objects, stlpath)

if __name__ == "__main__":
    myDC1 = synthetic()
    myDC1.createSynthetic(u"/home/nirb/regions/sin9.stl", 500, 500, 50, 50, 50, 50)

