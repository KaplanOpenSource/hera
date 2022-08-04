import unittest
import warnings
from hera_data_load_unittest import loader

class TestDatalayer(unittest.TestCase):
    def test_update_metadata_and_retrieve(self):
        warnings.simplefilter("ignore", ResourceWarning)
        source = "BNTL"
        dataType = "GIS_Buildings"
        projectName = 'testProject'
        expectedResult = f"Added source {source} to tool {dataType} in project {projectName}"
        foundResult  = f"Source {source} already exists in {projectName}"
        self.assertTrue(loader(projectName,dataType,source)==foundResult,expectedResult)

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
