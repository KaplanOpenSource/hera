import json
import os
import pandas
import numpy
from scipy.stats import norm
import  matplotlib.pyplot as plt 
from hera import datalayer
import unittest
import warnings

class TestAddingDataParquet(unittest.TestCase):
    
    def run_test_queries(self,y,projectName):
        # logging.getLogger(__name__).disabled = True

        x = numpy.linspace(norm.ppf(0.01), norm.ppf(0.99), 100)

        dataset1 = pandas.DataFrame(dict(x=x,y=norm.pdf(x,loc=0,scale=1)))
        dataset2 = pandas.DataFrame(dict(x=x,y=norm.pdf(x,loc=0,scale=0.5)))
        dataset3 = pandas.DataFrame(dict(x=x,y=norm.pdf(x,loc=0.5,scale=0.5)))

        workingdir = os.getcwd()
        # print(f"The current file directory is {workingdir}")

        dataset1File = os.path.join(workingdir,"dataset1.parquet")
        dataset2File = os.path.join(workingdir,"dataset2.parquet")
        dataset3File = os.path.join(workingdir,"dataset3.parquet")

        dataset1.to_parquet(dataset1File,engine='fastparquet',compression='GZIP')
        dataset2.to_parquet(dataset2File,engine='fastparquet',compression='GZIP')
        dataset3.to_parquet(dataset3File,engine='fastparquet',compression='GZIP')


        datalayer.Measurements.addDocument(projectName=projectName,
                                        type="Distribution",
                                        dataFormat=datalayer.datatypes.PARQUET,
                                        resource=dataset1File,
                                        desc=dict(loc=0,scale=1))

        datalayer.Measurements.addDocument(projectName=projectName,
                                        type="Distribution",
                                        dataFormat=datalayer.datatypes.PARQUET,
                                        resource=dataset2File,
                                        desc=dict(loc=0,scale=0.5))

        datalayer.Measurements.addDocument(projectName=projectName,
                                        type="Distribution",
                                        dataFormat=datalayer.datatypes.PARQUET,
                                        resource=dataset3File,
                                        desc=dict(loc=0.5,scale=0.5))
        
        # print(projectName)

        List1 = datalayer.Measurements.getDocuments(projectName=projectName,loc=y)
        return len(List1)


    def insert_retrieve(self,testID,projectName,x):
        warnings.simplefilter("ignore", ResourceWarning)

        dataset = pandas.DataFrame(dict(x=[x]))

        workingdir = os.getcwd()
        # print(f"The current file directory is {workingdir}")

        datasetFile = os.path.join(workingdir,f"{projectName}_{testID}.parquet")

        dataset.to_parquet(datasetFile,engine='fastparquet',compression='GZIP')
        datalayer.Measurements.addDocument(projectName=projectName,
                                type="Distribution",
                                dataFormat=datalayer.datatypes.PARQUET,
                                resource=datasetFile,
                                desc=dict(test_id = testID))

        result = datalayer.Measurements.getDocuments(projectName=projectName,test_id=testID)
        dataset_read = pandas.read_parquet(result[len(result)-1].resource)
        result_val = dataset_read.loc[0,'x']
        return result_val


    def insert_retrieve_dataframe(self,testID,projectName,dataset):
        warnings.simplefilter("ignore", ResourceWarning)

        # Get working direcotry
        workingdir = os.getcwd()
        datasetFile = os.path.join(workingdir,f"{projectName}_{testID}.parquet")

        # Exporting to parquet
        dataset.to_parquet(datasetFile,engine='fastparquet',compression='GZIP')
        datalayer.Measurements.addDocument(projectName=projectName,
                                type="Distribution",
                                dataFormat=datalayer.datatypes.PARQUET,
                                resource=datasetFile,
                                desc=dict(test_id = testID))

        # Retrieving data
        result = datalayer.Measurements.getDocuments(projectName=projectName,test_id=testID)
        dataset_read = pandas.read_parquet(result[len(result)-1].resource)
        return dataset_read

    # def test_query_count(self):
    #     # The exact assertion should return 2, for the nubmer of returnred queries
    #     returnedItems = self.run_test_queries(0,'test_1')
    #     self.assertEqual(returnedItems,2)
    

    def test_equal_returned_values(self):
        """
        Generates a random number and returns if the number from dataset is the same
        """
        number = numpy.random.randint(0,10000)
        result = self.insert_retrieve('test_exact_value','unittest_parquet',number)
        self.assertTrue(result==number,'Return number not equal to given number')



    def test_negative_nubmers(self):
        """
        The assertion should return True if the value is negative
        """
        number = numpy.random.randint(-1, -1000000)
        result = self.insert_retrieve('test_negative_values','unittest_parquet',number)
        self.assertLess(result,0,'Number not negative')

    def test_almost_equal(self):
        """
        The assertion should return a number close to the generated number
        """
        number = numpy.random.uniform(0.1, 1000000)
        number_close = number- 0.000001
        result = self.insert_retrieve('test_close_values','unittest_parquet',number_close)
        self.assertAlmostEqual(number, result,5,'Number not almost equal')

    def test_return_eqaul_lines(self):
        """
        The assertion should return a number rows of retrieved dataset, should be equal to inserted dataset number of rows
        """

        # Random number of rows
        r = numpy.random.randint(2,10)

        # Random List to populate rows
        random_list = numpy.random.random(r)

        df = pandas.DataFrame({'random_list':random_list})
        result = self.insert_retrieve_dataframe('test_eqaul_lines','unittest_parquet',df)
        self.assertEqual(r, len(result),'Number of rows not eqaul')
if __name__ =='__main__':
    unittest.main(warnings='ignore')
