import unittest
import hera_data_load

# continue here
#http://agiletesting.blogspot.com/2005/01/python-unit-testing-part-1-unittest.html

class TestDatalayer(unittest.TestCase):
    def test_update_metadata_and_retrieve(self):
        source = "BNTL"
        dataType = "GIS_Buildings"
        projectName = 'testProject'
        expectedResult = f"Added source {source} to tool {dataType} in project {projectName}"
        foundResult  = f"Source {source} already exists in {projectName}"
        self.assertTrue(hera_data_load.load(projectName,dataType,source)==foundResult,expectedResult)

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
