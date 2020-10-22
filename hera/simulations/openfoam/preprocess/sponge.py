import os
import pandas
import numpy
import json
from scipy.interpolate import interp1d
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from  PyFoam.Basics.DataStructures import Field
from ....utils import loggedObject

class spongeLayer(loggedObject):

    params=None

    _spongeLayersInterpolator = None

    def __init__(self,params,centerFile="0/C",scalarField="0/p"):

        super().__init__()

        self.logger.info("Initialize spongeLayer")

        self.logger.info("Loading configuration")
        if isinstance(params,str):
            self.logger.debug("params is a string, checking to see if it is a file on the disk")
            if os.path.exists(params):
                with open(params) as f:
                   self.params=json.load(f)
            else:
                self.params=json.loads(params)
        else:
            self.params=dict(params)

        self.logger.execution(f"got parameters {params}")
        self.logger.info("Reading cell centers")

        self._centerFile = ParsedParameterFile(centerFile)
        self._alpha      = ParsedParameterFile(scalarField)
        self._dirMap=dict(x=0,y=1,z=2)
        self._spongeLayers = {}

    def itrSpongeAreas(self):

        tmpField = Field([self.params['default_value']] * len(self._centerFile['internalField']), "alpha")

        for spongeName, spongeData in self.params["SpongeLayers"].keys():

            self.logger.execution(f"Creating sponge interpolation {spongeName}")
            coord = self._dirMap.get(spongeData['direction'])
            limits = pandas.DataFrame.from_dict(spongeData['vals']).sort_values(self.params['params']['direction'])

            tmp=calcSpongeArea(tmpField,coord,limits)
            #_spongeLayers['f = interp1d(limits[coord], limits['alpha'], kind=type)

    def calcSpongeArea(self,tmpField,coord,limits):
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
        self.logger.info(f"Starting alpha Mat with type {type}")
        f = interp1d(limits[coord], limits['alpha'], kind=type)
        for index,point in enumerate() self._centerFile['internalField']):
            if index%50000==0:
                self.logger.execution(f"processing node {index} of {len(self._centerFile['internalField'])}")
                self.logger.execution(f"Processing sponge {spongeName}")

                val = point[coord]
                try:
                    tmpField[index] = self.get_alpha(f, val)
                except:
                    'log here'




    def CreateAlphaFile(self,f,val):

        tmpField=itrSpongeAreas()
        self._alpha.content['dimensions']= self.params['dimensions']
        self._alpha['internalField']=tmpField

        self._alpha.writeFileAs(os.path.join(self.params['fileName']))


