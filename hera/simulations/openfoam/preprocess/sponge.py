import os
import pandas
import numpy
import json
from scipy.interpolate import interp1d
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from  PyFoam.Basics.DataStructures import Field
from ....utils import loggedObject

class spongeLayer(loggedObject):
    """
    This class manages the creation of an alpha file which will be used in the run to create a sponge layer.
    The class receives 'params' as a dict, string or path to json file with the parameters in the following structure:

    filename- the filename to save on disc
    default_value- the initial fill value to the temporary internalField which will last in the area out of any sponge area
    dimensions- the dimensions of the alpha parameter
    SpongeLayers- the sponge layers to aplly
        "sponeName":"profile" : "TypeOfProfile", "direction": "x or y or z","vals": [{"direction string": diirection value,"alpha": alpha _value}


    params example:

    {
  "fileName": "alpha",
  "default_value" : 0,
  "dimensions": "[ 0 0 -1 0 0 0 0 ]",
  "SpongeLayers" : {
    "top": {
      "profile": "linear",
      "direction": "z",
      "vals":
      [{"z": 7,"alpha": 0},
       {"z": 8,"alpha": 1e-2},
       {"z": 9,"alpha": 1e-1},
       {"z": 10,"alpha": 1e-1}]
        }
      }
    }


    """

    params      = None
    _centerFile = None
    _alpha      = None
    _dirMap     = None

    def __init__(self, params, CellCentersFile="0/C", ScalarFieldFile="0/p"):

        """

        parametrs
        --------
        params: The parameters
        CellCentersFile: the path to cell center file
        ScalarFieldFile: the path to scalar tamplate file
        """

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

        self.logger.execution(f"got parameters from {params}")
        self.logger.info("Reading cell centers")

        self._centerFile = ParsedParameterFile(CellCentersFile)
        self.logger.info("Reading alpha template")
        self._alpha      = ParsedParameterFile(ScalarFieldFile)
        self._dirMap=dict(x=0,y=1,z=2)


    def itrSpongeAreas(self):

        """
        iteration on sponge areas:
            - Creates a temporary 'internalField'
            - looping on all sponge areas from params file

        :return:
        """

        tmpField = Field([self.params['default_value']] * len(self._centerFile['internalField']), "List<scalar>")  #"alpha"

        for spongeName, spongeData in self.params["SpongeLayers"].items():

            self.logger.execution(f"Creating sponge interpolation {spongeName}")

            coord = spongeData['direction']
            coord_val= self._dirMap.get(spongeData['direction'])
            limits = pandas.DataFrame.from_dict(spongeData['vals']).sort_values(spongeData['direction'])
            type=spongeData['profile']

            tmp=self.calcSpongeArea(tmpField,coord,coord_val,limits,spongeName,type)

            return tmp

    def calcSpongeArea(self,tmpField,coord,coord_val,limits,spongeName,type):

        """
        fills the tmp matrix of internalField with alpha values according to the params

        parameters
        ----------

        tmpField:
        coord:
        coord_val:
        limits:
        spongeName:
        type: the interpolation method (from params file)

        return
        ------

        tmpField: the new alpha internalField
        """


        self.logger.info(f"Starting alpha Mat with type {type}")

        f = interp1d(limits[coord], limits['alpha'], kind=type)
        numOfCells=len(self._centerFile['internalField'])

        self.logger.execution(f"Processing sponge {spongeName}")
        for index,point in enumerate(self._centerFile['internalField']):
            if index / numOfCells * 100 % 10 == 0:

                self.logger.execution(f"processing node {index} of {len(self._centerFile['internalField'])}")

            val = point[coord_val]
            try:
                tmpField[index] = f(val)
            except:
                self.logger.debug(f"point {point} not in sponge {spongeName}")
        self.logger.execution(f"Processing sponge {spongeName} finished")

        return tmpField


    def CreateAlphaFile(self):

        """
        Main flow of Alpha file creation.
        Runs all necessary functions from the spongeLayer class and saves the result as a 'fileName' file

        """

        tmpField = self.itrSpongeAreas()

        self._alpha.content['dimensions']= self.params['dimensions']
        self._alpha['internalField']=tmpField

        self._alpha.writeFileAs(os.path.join('constant',self.params['fileName']))


