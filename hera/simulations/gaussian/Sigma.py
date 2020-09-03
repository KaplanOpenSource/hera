import pandas
import numpy
from ..utils import toUnum,toNumber
from unum.units import *


class BriggsRural(object):

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

    def getSigma(self,x,stability):

        x = numpy.array([toUnum(y,m).asNumber() for y in numpy.atleast_1d(x)])
        Ax, Bx, Cx = self._coeffX.loc[stability][['A','B','C']]
        Az, Bz, Cz = self._coeffZ.loc[stability][['A', 'B', 'C']]

        return pandas.DataFrame({
                'sigmaX' : Ax*x*(1+Bx*x)**Cx,
                'sigmaZ' : Az*x*(1+Bz*x)**Cz,'distance' : x})





briggsRural = BriggsRural()
