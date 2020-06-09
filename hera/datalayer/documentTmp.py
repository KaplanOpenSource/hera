import pandas

class AbstractDocument(object):
    _data = None

    def __init__(self, data):
        self._data = data

    @property
    def type(self):
        return self._data['type']

    @property
    def name(self):
        return self._data['name']



class GIS_Document(AbstractDocument):
    def __init__(self, data):
        super.__init__(data=data)



class Experimental_Document(AbstractDocument):
    def __init__(self, data):
        super.__init__(data=data)
    
    def getStations(self):
        """

        :return: Stations pandas

        TODO:
            - Change to a general desc.

        """
        return pandas.DataFrame(self._data['stations'])

    def getDevices(self):
        """

        :return: Devices pandas
        """
        return pandas.DataFrame(self._data['devices'])


class Numerical_Document(AbstractDocument):
    def __init__(self, data):
        super.__init__(data=data)

class Analysis_Document(AbstractDocument):
    def __init__(self, data):
        super.__init__(data=data)
