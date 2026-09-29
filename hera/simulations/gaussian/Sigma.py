import pandas
import numpy
from hera.utils import tounit,tonumber
from hera.utils.unitHandler import ureg, unumToPint
from scipy import optimize


class AbstractSigma:

    def getVirtualDistance(self,sigma0,stability):
        """
        Calculates the virtual distances for a given sigma0.

        Parameters
        ----------
        sigma0 : 3-tuple of float/unum (default m)
            The initial cloud size (standard deviation) in the x,y and z dimensions.

        stability : str
            Must be A-F (capital letters)

        Returns
        -------

        """

        # 1. use root finding to find Ix and use self.getSigma(x,stability)
        # 2. use root finding to find Iy and use self.getSigma(x,stability)
        # 3. use root finding to find Iz and use self.getSigma(x,stability)

        stability=stability

        sx,sy,sz = sigma0
        sx = tonumber(sx, ureg.m)
        sy = tonumber(sy, ureg.m)
        sz = tonumber(sz, ureg.m)

        sigma_x = lambda x: self.getSigma(x=x, stability=stability, units=False)['sigmaX'][0] - sx
        sigma_y = lambda x: self.getSigma(x=x, stability=stability, units=False)['sigmaY'][0] - sy
        sigma_z = lambda x: self.getSigma(x=x, stability=stability, units=False)['sigmaZ'][0] - sz

        Ix = optimize.newton(sigma_x, x0=0)*ureg.m
        Iy = optimize.newton(sigma_y, x0=0) * ureg.m
        Iz = optimize.newton(sigma_z, x0=0) * ureg.m

        return Ix, Iy, Iz

    def getSigma(self,x,stability,sigma0=None, **kwargs):
        """
            Return the sigma.

            Implemented in the right sigma class.
        Parameters
        ----------
        x
        stability
        sigma0

        Returns
        -------

        """
        raise NotImplementedError("The function getSigma is not implemented for abstractSigma. Use the son. ")


class BriggsRural(AbstractSigma):

    _coeffX = None
    _coeffZ = None

    def __init__(self):

        self._coeffX = pandas.DataFrame({
            'A' : [0.22,0.16,0.11,0.08,0.06,0.04],
            'B' : [1e-4]*6,
            'C' : [-0.5]*6},
            index=['A','B','C','D','E','F'])

        self._coeffZ = pandas.DataFrame({
            'A' : [0.2,0.12,0.08,0.06,0.03,0.016],
            'B' : [0,0,2e-4,1.5e-3,3e-4,3e-4],
            'C' : [1,1,-0.5,-0.5,-1,-1]},
            index=['A','B','C','D','E','F'])


    def __call__(self,x,stability):
        return self.getSigma(x,stability)

    def getSigma(self, x, stability, sigma0=None, units=True):
        """
            Computes the briggs sigma and return the sigma for the request
            distances in the stability. Taking the initial size to be sigma0.

        Parameters
        ----------
        x : numpy array/list/float/ unum [default units m]
            The distance from the source.

        stability : str
            Must be A-F (capital letters)

        sigma0 : 3-tuple of float/unum (default m)
            The initial cloud size in the x,y and z dimensions.
            If None, use a point source.

        Returns
        -------
            pandas.DataFrame


        """

        x = numpy.array([tonumber(y,ureg.m) for y in numpy.atleast_1d(x)])
        Ax, Bx, Cx = self._coeffX.loc[stability][['A','B','C']]
        Az, Bz, Cz = self._coeffZ.loc[stability][['A', 'B', 'C']]

        # Compute Ix,Iy,Iz
        if sigma0 is None:
            Ix = Iy = Iz = 0
        else:
            Ix,Iy,Iz = self.getVirtualDistance(sigma0,stability)


        Ix = tonumber(Ix,ureg.m)
        Iy = tonumber(Iy,ureg.m)
        Iz = tonumber(Iz,ureg.m)

        if units:
            dict_res = {
                'sigmaX': Ax * (x + Ix) * (1 + Bx * (x + Ix)) ** Cx * ureg.m,
                'sigmaY': Ax * (x + Iy) * (1 + Bx * (x + Iy)) ** Cx * ureg.m,
                'sigmaZ': Az * (x + Iz) * (1 + Bz * (x + Iz)) ** Cz * ureg.m, 'distance': x * ureg.m
            }

        else:
            dict_res = {
                'sigmaX': Ax * (x + Ix) * (1 + Bx * (x + Ix)) ** Cx,
                'sigmaY': Ax * (x + Iy) * (1 + Bx * (x + Iy)) ** Cx,
                'sigmaZ': Az * (x + Iz) * (1 + Bz * (x + Iz)) ** Cz, 'distance': x
            }


        return dict_res




briggsRural = BriggsRural()
