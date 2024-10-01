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

    def test_01_addRepositories(self):
        self._addRepositories()
        print("\nutils.data.toolkit.dataToolkit.addRepository function tested successfully.\n")

    def test_02_getRepository(self):
        test_cases = self.json_file['test_addRepository']
        self._addRepositories()
        for case in test_cases:
            with self.subTest(case=case):
                rep = self.tk.getRepository(case['repositoryName'])
                self.assertTrue(rep,f"Could not use getRepository to Repository:{case['repositoryName']}")

        print("\nutils.data.toolkit.dataToolkit.getRepository function tested successfully.\n")

    def test_03_loadAllDatasourcesInRepositoryToProject(self):
        test_cases = self.json_file['test_loadAllDatasourcesInRepositoryToProject']
        self._addRepositories()
        for case in test_cases:
            with self.subTest(case=case):
                proj = Project(case['projectName'])
                self.tk.loadAllDatasourcesInRepositoryToProject(proj.projectName,case['repositoryName'],overwrite=True)
                count_loaded_to_project = len([doc.to_mongo().to_dict() for doc in proj.getMeasurementsDocuments()])

                count = 0
                with open(os.path.join(os.getcwd(), "testCases", case['repositoryName']+".json")) as file:
                    f = json.load(file)

                for key, value in f.items():
                    if 'DataSource' in value:
                        count += len(value['DataSource'])
                    if 'Measurements' in value:
                        count += len(value['Measurements'])

                self.assertEqual(count_loaded_to_project,count,f"Number of Repository documents {count} is not equal to projects documents size {count_loaded_to_project}")

        print("\nutils.data.toolkit.dataToolkit.loadAllDatasourcesInRepositoryToProject function tested successfully.\n")

    def _addRepositories(self):
        test_cases = self.json_file['test_addRepository']
        for case in test_cases:
            repository_path = os.path.join(os.getcwd(), "testCases", case['repositoryPath'])
            self.tk.addRepository(repositoryName=case['repositoryName'], repositoryPath=repository_path,
                                  overwrite=True)
    @classmethod
    def tearDown(self):
        for project_name in self.json_file["projects_to_delete"]:
            proj = Project(project_name)
            [x.delete() for x in proj.getMeasurementsDocuments()]

        for rep in self.json_file["test_addRepository"]:
            self.tk.deleteDataSource(datasourceName=rep['repositoryName'])


if __name__ == '__main__':
    unittest.main()
