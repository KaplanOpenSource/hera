import unittest
import json
from hera.utils.data.toolkit import dataToolkit
from hera.datalayer.project import Project
import os
class TestRepository(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        with open('testCases/repository.json', 'r') as file:
            self.json_file = json.load(file)

        self.tk = dataToolkit()

    def test_01_addRepository(self):
        test_cases = self.json_file['test_addRepository']
        for case in test_cases:
            with self.subTest(case=case):
                repository_path = os.path.join(os.getcwd(),"testCases",case['repositoryPath'])
                self.tk.addRepository(repositoryName=case['repositoryName'],repositoryPath=repository_path,overwrite=True)

        print("\nutils.data.toolkit.dataToolkit.addRepository function tested successfully.\n")

    def test_02_getRepository(self):
        test_cases = self.json_file['test_addRepository']
        for case in test_cases:
            with self.subTest(case=case):
                rep = self.tk.getRepository(case['repositoryName'])
                self.assertTrue(rep,f"Could not use getRepository to Repository:{case['repositoryName']}")

        print("\nutils.data.toolkit.dataToolkit.getRepository function tested successfully.\n")

    def test_03_loadAllDatasourcesInRepositoryToProject(self):
        test_cases = self.json_file['test_addRepository']
        for case in test_cases:
            with self.subTest(case=case):
                proj = Project(case['projectName'])
                self.tk.loadAllDatasourcesInRepositoryToProject(proj.projectName,case['repositoryName'],overwrite=True)





if __name__ == '__main__':
    unittest.main()
