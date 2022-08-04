import unittest
from hera_data_load_unittest import loader
import json
import pathlib
import os
import warnings

class Test_Datalayer(unittest.TestCase):
    def setUp(self):
        f = open(os.path.join(pathlib.Path(__file__).parent,"datasources.json"))
        self.json_file = json.load(f)
        f.close()
    def test_update_metadata_and_retrieve(self):
        warnings.simplefilter("ignore", ResourceWarning)
        projectName = 'Test Project 1'
        for dataType in self.json_file['dataType'].keys():
            for source in self.json_file['dataType'][dataType].keys():
                with self.subTest(f'----------{projectName} {dataType} {source}---------'):
                    expectedResult = f"Added source {source} to tool {dataType} in project {projectName}"
                    foundResult  = f"Source {source} already exists in {projectName}"
                    self.assertTrue(loader(projectName,dataType,source)==foundResult,expectedResult)   

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'],warnings='ignore')
