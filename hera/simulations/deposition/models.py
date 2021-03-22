from ...datalayer import project
from ...riskassessment import AgentHome
from ..utils import toNumber, toUnum
from unum.units import *
import numpy

class depositionModels(object):

    _surface = None
    _depositionModel = None
    _ustar = None
    _density = None
    _temperature = None
    _diameter = None
    _heatFlux = None

    @property
    def surface(self):
        return self._surface

    @surface.setter
    def surface(self,newSurface):
        self._surface=newSurface

    @property
    def ustar(self):
        return self._ustar

    @ustar.setter
    def ustar(self,newustar):
        self._ustar=newustar

    @property
    def heatFlux(self):
        return self._heatFlux

    @heatFlux.setter
    def heatFlux(self,newheatFlux):
        self._ustar=newheatFlux

    @property
    def diameter(self):
        return self._diameter

    @diameter.setter
    def ustar(self,newdiameter):
        self._diameter=newdiameter

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self,newtemperature):
        self._temperature=newtemperature

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self,newdensity):
        self._density=newdensity

    @property
    def depositionModel(self):
        return self._depositionModel

    @depositionModel.setter
    def depositionModel(self,newModel):
        self._depositionModel=newModel

    def __init__(self, projectName = "deposition", surface={"name":"Desert","zrough":0.04},
                 ustar=1, density=1500, diameter = 1E-6*m,heatFlux=0.1*W/(m**2), temperature=293*K, depositionModel="Petroff"):

        p = project.Project(projectName)
        self._depositionModel = depositionModel
        self._ustar = ustar
        self._density = density
        self._temperature = toNumber(toUnum(temperature,K),K)
        self._surface = surface if type(surface)==dict else p.getCacheDocuments(type="surface",surface=surface)[0].asDict()["desc"]
        self._diameter = toNumber(toUnum(diameter,m),m)
        self._heatFlux = toNumber(toUnum(heatFlux, W/(m**2)), W/(m**2))

    def depositionRate(self):
        return getattr(self, f"depositionRate_{self._depositionModel}")()

    def depositionRate_Petroff(self):

        nuA = 0.0000157
        muA = 0.0000189
        g = 9.81
        lpm = 0.000000067
        kb = 1.83E-23
        kappa = 0.4
        rhoA = 1.2
        cpA = 1000

        ustar = self.ustar
        rhop = self.density
        Tpe = self.temperature
        dp = self.diameter
        H = self.heatFlux

        zrough = self.surface["zrough"] if "zrough" in self.surface.keys() else 0  # roughness length
        hcanop = self.surface["hcanop"] if "hcanop" in self.surface.keys() else 0  # height of canopy
        hdepl = self.surface["displacementHeight"] if "displacementHeight" in self.surface.keys() else 0
        obstacleShape = self.surface["obstacleShape"] if "obstacleShape" in self.surface.keys() else None

        def PhiH(xi):
            if xi >= -2 and xi < 0:
                ret = (1 - 16 * xi) ** (-0.5)
            elif xi >= 0 and xi <= 1:
                ret = 1 + 5 * xi
            else:
                raise KeyError("PhiH is out of stability range!")
            return ret

        dpm = 0.000001*dp
        cu = 1+2*lpm/dpm*(1.257+0.4*numpy.exp(-1.1*dpm/(2*lpm))) # Cunningham slip factor
        Db = (cu*kb*Tpe)/(3*numpy.pi*muA*dpm)# Brownian diffusivity
        Sc = nuA/Db # Schmidt number
        trel = rhop*dpm**2*cu/(18*muA) # Particle relaxation time
        egb = Sc ** (-2 / 3) / (14.5 * (numpy.pi / (6 * numpy.sqrt(3)) + 1 / numpy.sqrt(3) * numpy.arctan(
              (2 * Sc ** (1 / 3) / 2.9 - 1) / numpy.sqrt(3)) +1 / 6 * numpy.log((1 + (Sc ** (1 / 3) / 2.9)) ** 2 /
              (1 - Sc ** (1 / 3) / 2.9 + (Sc ** (1 / 3) / 2.9) ** 2))))  # ground deposition
        LO = -ustar ** 3 / kappa * Tpe / g * rhoA * cpA / H
        lm = kappa * (hcanop - hdepl) / PhiH((hcanop - hdepl) / LO)

        if obstacleShape is not None:

            def PhiM(xi):
                if xi >= -2 and xi < 0:
                    ret = (1 - 16 * xi) ** (-0.25)
                elif xi >= 0 and xi <= 1:
                    ret = 1 + 5 * xi
                else:
                    raise KeyError("PhiM is out of stability range!")
                return ret

            def PsiM(xi):
                if xi >= -2 and xi < 0:
                    ret = 2*numpy.log(0.5*(1+(1 - 16 * xi) ** (0.25)))+numpy.log(0.5*(1+(1 - 16 * xi) ** (0.5)))
                elif xi >= 0 and xi <= 1:
                    ret = -5 * xi
                else:
                    raise KeyError("PsiM is out of stability range!")
                return ret

            lai = self.surface["LAI"]  # Leaf area index
            kx = self.surface["plagiophile"]
            xb = self.surface["xb"]
            xim = self.surface["xim"]
            xin = self.surface["xin"]
            xit = self.surface["xit"]
            ObstacleSize = self.surface["ObstacleSize"]

            Uh = ustar / kappa * numpy.log((hcanop - hdepl) / zrough) - PsiM((hcanop - hdepl) / LO) + PsiM(zrough / LO)
            alphaVeg = (kx * lai / (12 * kappa ** 2 * (1 - hdepl / hcanop) ** 2)) ** (1 / 3) * (PhiM((hcanop - hdepl) / LO)) ** (2 / 3)
            egt = (0.00035 / (nuA ** 2)) * trel ** 2 * (ustar * numpy.exp(-alphaVeg)) ** 4 if trel * (
                    ustar * numpy.exp(-alphaVeg)) ** 2 / nuA < 20 else 0.14

            eb = {"leaf": (xb * 0.664 * (nuA / Db) ** (-2 / 3) * (Uh * ObstacleSize / nuA) ** (-0.5)),
                  "needle":(xb * 0.467 * (nuA / Db) ** (-2 / 3) * (Uh * ObstacleSize / nuA) ** (-0.5))}
            ein = {"leaf": (xin*kx/2*dp*0.000001/ObstacleSize*(2+numpy.log(4*ObstacleSize/(dp*0.000001))))}
            eim = {"leaf": (xim*kx/(1+0.47*ObstacleSize/(trel*Uh))**2)}
            eit = {"leaf": xit*0.14}
            Qveg = hcanop/lm*lai*(Uh/ustar*eb[obstacleShape]+ein[obstacleShape]+eim[obstacleShape]+eit[obstacleShape])
            Qsol = hcanop/lm*(egb+egt)
            eta = numpy.sqrt(1+4*Qveg/alphaVeg**2)

            vds = ustar*(egb+egt)*(1+(Qveg/Qsol-alphaVeg/2)*2/(alphaVeg*eta)*numpy.tanh(alphaVeg*eta/2))/(1+(Qsol+alphaVeg/2)*2/(alphaVeg*eta)*numpy.tanh(alphaVeg*eta/2)) # deposition rate
        else:
            Qsol = hcanop / lm *egb if lm != 0 else 0
            vds = ustar*egb/(1+Qsol)
        return vds
