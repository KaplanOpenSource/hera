import numpy
import pandas

class sourcesFactoryTool():

    @property
    def sourcesTypeList(self):
        return [x.split("_")[1] for x in dir(self) if "makeSource" in x and "_" in x]

    def makeSource(self, x, y, z, nParticles, type="Point", **kwargs):

        slist = self.sourcesTypeList
        if type not in slist:
            raise ValueError(f"The type must be [{','.join(slist)}]. Got {type} instead")

        return getattr(self, f"makeSource_{type}")(x, y, z, nParticles, **kwargs)

    def makeSource_Point(self,x,y,z,nParticles,**kwargs):
        return pandas.DataFrame({"x":[x for i in range(nParticles)],"y":[y for i in range(nParticles)],"z":[z for i in range(nParticles)]})

    def makeSource_Circle(self,x,y,z,nParticles,radius,distribution="uniform",**kwargs):

        Rs = getattr(numpy.random,distribution)(0,radius,nParticles)
        thetas = numpy.random.uniform(0,2*numpy.pi,nParticles)
        xs = []
        ys = []
        for R, theta in zip(Rs,thetas):
            xs.append(x+R*numpy.cos(theta))
            ys.append(y+R*numpy.sin(theta))
        return pandas.DataFrame({"x":xs,"y":ys,"z":[z for i in range(nParticles)]})

    def makeSource_Sphere(self,x,y,z,nParticles,radius,distribution="uniform",**kwargs):

        Rs = getattr(numpy.random,distribution)(0,radius,nParticles)
        thetas = numpy.random.uniform(0,2*numpy.pi,nParticles)
        phis = numpy.random.uniform(0,2*numpy.pi, nParticles)
        xs = []
        ys = []
        zs = []
        for R, theta, phi in zip(Rs,thetas,phis):
            xs.append(x+R*numpy.sin(theta)*numpy.cos(phi))
            ys.append(y+R*numpy.sin(theta)*numpy.sin(phi))
            zs.append(z+R*numpy.cos(theta))
        return pandas.DataFrame({"x":xs,"y":ys,"z":zs})

    def makeSource_Cylinder(self,x,y,z,nParticles,radius,height,horizontalDistribution="uniform",verticalDistribution="uniform",**kwargs):

        Rs = getattr(numpy.random,horizontalDistribution)(0,radius,nParticles)
        Hs = getattr(numpy.random,verticalDistribution)(-height/2,height/2,nParticles)
        thetas = numpy.random.uniform(0,2*numpy.pi,nParticles)
        xs = []
        ys = []
        zs = []
        for R, theta,H in zip(Rs,thetas,Hs):
            xs.append(x+R*numpy.cos(theta))
            ys.append(y+R*numpy.sin(theta))
            zs.append(z+H)
        return pandas.DataFrame({"x":xs,"y":ys,"z":zs})

    def makeSource_Rectangle(self,x,y,z,nParticles,lengthX,lengthY,rotateAngle=0,**kwargs):

        xdist = numpy.random.uniform(-lengthX/2,lengthX/2,nParticles)
        ydist = numpy.random.uniform(-lengthY/2, lengthY/2, nParticles)
        xs = []
        ys = []
        for i in range(nParticles):
            xs.append(x+xdist[i]*numpy.cos(rotateAngle)+ydist[i]*numpy.sin(rotateAngle))
            ys.append(y-xdist[i]*numpy.sin(rotateAngle)+ydist[i]*numpy.cos(rotateAngle))
        return pandas.DataFrame({"x":xs,"y":ys,"z":[z for i in range(nParticles)]})

    def makeSource_Cube(self,x,y,z,nParticles,lengthX,lengthY,lengthZ,rotateAngle=0,**kwargs):

        xdist = numpy.random.uniform(-lengthX/2,lengthX/2,nParticles)
        ydist = numpy.random.uniform(-lengthY/2,lengthY/2, nParticles)
        zdist = numpy.random.uniform(-lengthZ/2,lengthZ/2, nParticles)
        xs = []
        ys = []
        zs = []
        for i in range(nParticles):
            xs.append(x+xdist[i]*numpy.cos(rotateAngle)+ydist[i]*numpy.sin(rotateAngle))
            ys.append(y-xdist[i]*numpy.sin(rotateAngle)+ydist[i]*numpy.cos(rotateAngle))
            zs.append(z+zdist[i])
        return pandas.DataFrame({"x":xs,"y":ys,"z":zs})