import unittest
import hera_data_load
import json
import pathlib
import os

class Test_Datalayer(unittest.TestCase):
    def setUp(self):
        f = open(os.path.join(pathlib.Path(__file__).parent,"datasources.json"))
        self.json_file = json.load(f)
        f.close()
    def test_update_metadata_and_retrieve(self):
        projectName = 'Test Project'
        for dataType in self.json_file['dataType'].keys():
            for source in self.json_file['dataType'][dataType].keys():
                with self.subTest(f'{projectName} {dataType} {source}'):
                    expectedResult = f"Added source {source} to tool {dataType} in project {projectName}"
                    foundResult  = f"Source {source} already exists in {projectName}"
                    self.assertTrue(hera_data_load.load(projectName,dataType,source)==foundResult,expectedResult)
    

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
