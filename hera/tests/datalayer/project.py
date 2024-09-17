import unittest
import json
from hera.datalayer.project import createProjectDirectory
from hera.datalayer.project import Project
class TestProject(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        with open('testCases/project.json', 'r') as file:
            self.json_file = json.load(file)

        self.projects = {}

        test_cases = self.json_file['__init__']
        for case in test_cases:
            proj = Project(case['projectName'])
            self.projects[case['projectName']] = proj
        print(f"\nsetUpClass: Projects initialized: {self.projects.keys()}\n")

    def test_01_init(self):
        for project_name, proj in self.projects.items():
            self.assertTrue(proj, f"Project {project_name} initialization failed.")
        print("\nhera.datalayer.project.Project.__init()__ function tested successfully.\n")

    def test_02_addMeasurementsDocument(self):
        test_cases = self.json_file['test_addMeasurementsDocument']
        for case in test_cases:
            with self.subTest(case=case):
                proj = self.projects[case['projectName']]
                doc = proj.addMeasurementsDocument(resource=case['resource'],
                                                                           dataFormat=case['dataFormat'],
                                                                           type=case['type'],
                                                                           desc=case['desc'])
                self.assertTrue(doc,f"Document {case} was not added to project {case['projectName']}")

        print("\ndatalayer.project.Project.addMeasurementsDocument function tested successfully.\n")

    def test_03_getMeasurementsDocuments(self):
        test_cases = self.json_file['test_addMeasurementsDocument']
        for case in test_cases:
            with self.subTest(case=case):
                proj = self.projects[case['projectName']]
                doc = proj.getMeasurementsDocuments(resource=case['resource'],
                                                   dataFormat=case['dataFormat'],
                                                   type=case['type'],
                                                   *case['desc'])
                self.assertNotEqual(list(doc),[],f"getMeasurementsDocuments in {case} was not performed correctly.")

        print("\ndatalayer.project.Project.getMeasurementsDocuments function tested successfully.\n")

    def test_04_deleteMeasurementsDocuments(self):
        test_cases = self.json_file['test_deleteMeasurementsDocuments']
        for case in test_cases:
            with self.subTest(case=case):
                proj = self.projects[case['projectName']]
                doc = proj.deleteMeasurementsDocuments()
                self.assertNotEqual(list(doc), [], f"getMeasurementsDocuments in {case} was not performed correctly.")

    print("\ndatalayer.project.Project.deleteMeasurementsDocuments function tested successfully.\n")

    @classmethod
    def tearDownClass(self):
        for proj in self.projects.values():
            [x.delete() for x in proj.getMeasurementsDocuments()]

if __name__ == '__main__':
    unittest.main()