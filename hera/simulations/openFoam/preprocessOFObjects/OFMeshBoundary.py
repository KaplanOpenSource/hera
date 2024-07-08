
#########################################################################
#
#               Mesh handling
#########################################################################
class OFMeshBoundary:
    _case = None
    _checkIfParallel = None

    _boundaryNames = None

    @property
    def case(self):
        return self._case

    def __init__(self, directory: str, checkParallel: bool = True):
        """
            Reads the boundary from the boundary file of the mesh.

            If checkParallel is true, then checks if the parallel case exists. If it does,
            read the boundary from the parallel case. Else, read it from the constant.



        Parameters
        ----------
        baseFile
        """
        self._case = directory
        self._checkIfParallel = checkParallel

        self._boundaryNames = []
        if self._checkIfParallel and os.path.exists(os.path.join(self.case, "processor0")):
            for proc in glob.glob(os.path.join(self.case, "processor*")):
                self._boundaryNames += self._readBoundary(os.path.join(proc, "constant", "polyMesh", "boundary"))
        else:
            self._boundaryNames = self._readBoundary(os.path.join(directory, "constant", "poly", "boundary"))

    def getBoundary(self, filterProcessor: bool = True):
        """
            Return a list of all the boundaries. If fileterprocessor is true,
            remove all the processor*.
        Parameters
        ----------
        filterProcessor : bool
            If true remove all the processor* faces from the list.

        Returns
        -------

        """
        return list(set([x for x in self._boundaryNames if
                         'procBoundary' not in x] if filterProcessor else self._boundaryNames))

    def _readBoundary(self, boundaryFile):
        """
                Reads the boundary file and extracts the boundary names.
        Parameters
        ----------
        boundaryFile

        Returns
        -------

        """

        def isInt(line):
            try:
                return int(line)
            except:
                return None

        with open(boundaryFile, "r") as inFile:
            data = inFile.readlines()

        firstLine = [i for i, x in enumerate(data) if isInt(x) is not None][0]
        braceList = [i for i, x in enumerate(data[firstLine:]) if x.strip() == '{']
        return [data[firstLine:][x - 1].strip() for x in braceList]

