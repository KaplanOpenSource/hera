import os
import pandas
import numpy
import json
from scipy.interpolate import interp1d
from hera.simulations.openfoam.utils import centersToPandas as CtP


class spongeLayer(object):

    params=None

    def __init__(self,params):

        name =  ".".join(str(self.__class__)[8:-2].split(".")[1:])
        self._logger = logging.getLogger(name)


        if isinstance(params,str):
            if os.path.exists(params):

                with open(params) as f:
                   self.params=json.load(f)
            else:
                self.params=json.loads(params)
        else:
            self.params=dict(params)
        print(self.params)

    def alphaMat(self,type='linear'):
        """
        builds matrix of [x,y,z,alpha]

        parameters
        ----------
        type: the unterpolation method

        return
        ------
        P: dataframe
            the dataframe with the retrived alpha and coordinates
        """

        Cpath=self.params["Cpath"]
        coor = self.params['params']['direction']

        # read the cell centers data
        P=CtP(skipend=self.params['skipend'],
              filepath=Cpath,
              skiphead=self.params['skiphead'])

        P['alpha']=numpy.nan

        limits=pandas.DataFrame.from_dict(self.params['params']['vals'])\
                               .sort_values(self.params['params']['direction'])

        f = interp1d(limits[coor], limits['alpha'], kind=type)

        for index,data in P.iterrows():
            if index%50000==0:
                print(index)
            val=data[self.params['params']['direction']]
            P['alpha'][index]=self.get_alpha(f,val)

        return P


    # def get_alpha(self,limits,val,coor='z',type='linear'):
    #     """
    #     parameters
    #     ----------
    #     limits: dataframe
    #         the values to interpolate between
    #     val: float
    #         the specific value in 'direction' to interpolate
    #     coor: the direction of interpolation
    #     type: the method of interpolation
    #     return
    #     ------
    #     alpha: the resulted value
    #
    #     """
    #
    #     if val<= self.params['params']['start']:
    #         alpha=self.params['params']['default_value']
    #     else:
    #         lower_coor=limits[coor][limits[coor]<val].max()
    #         upper_coor=limits[coor][limits[coor]>val].min()
    #         ar = [float(limits['alpha'].loc[limits[coor]==lower_coor]),
    #               float(limits['alpha'].loc[limits[coor]==upper_coor])]
    #         f = interp1d([lower_coor,upper_coor], ar, kind=type)
    #         alpha=f(val)
    #
    #     return alpha



    def get_alpha(self,f,val):
        """
        parameters
        ----------
        f: dataframe
            the interpolation function
        val: float
            the specific value in 'direction' to interpolate

        return
        ------
        alpha: the resulted value

        """

        if val<= self.params['params']['start']:
            alpha=self.params['params']['default_value']
        else:
            alpha=f(val)

        return alpha

    def makeAlphaFile(self):


        header= """
        /*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      scalarField;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];
        """

        pass