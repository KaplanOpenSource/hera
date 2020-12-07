import logging


class loggedObject:

    _logger = None

    @property
    def logger(self):
        return self._logger

    def __init__(self):
        name = ".".join(str(self.__class__)[8:-2].split(".")[1:])
        self._logger = logging.getLogger(name)