import os
import pandas
import numpy
from scipy.stats import norm
from hera import datalayer
import unittest
import warnings

class TestDatalayerParquet(unittest.TestCase):
    
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
        dataset_read = pandas.read_parquet(result[0].resource)
        return dataset_read

    def overwrite_parquet(self,testID,projectName,newDataset):
        """
        function retreives location of dataset via datalayer and overwrites existing value, with a new dataset and returns the new read parquet file
        """
        warnings.simplefilter("ignore", ResourceWarning)
        to_update = datalayer.Measurements.getDocuments(projectName=projectName,test_id=testID)
        datasetFille =to_update[0].resource
        newDataset.to_parquet(datasetFille,engine='fastparquet',compression='GZIP')
        
        # Retrieving data
        dataset_read = pandas.read_parquet(datasetFille)
        return dataset_read

    # def test_query_count(self):
    #     # The exact assertion should return 2, for the nubmer of returnred queries
    #     returnedItems = self.run_test_queries(0,'test_1')
    #     self.assertEqual(returnedItems,2)
    
    def setUp(self):
        print('setUp')
        datalayer.Measurements.deleteDocuments(projectName='unittest_parquet')     

    def tearDown(self):
        print('tearDown')
        datalayer.Measurements.deleteDocuments(projectName='unittest_parquet')
        # Get working direcotry
        workingdir = os.getcwd()
        files = os.listdir(workingdir)
        files = [file for file in files if ".parquet" in file]
        for file in files:
            if file.split('.')[1]=='parquet':
                os.remove(os.path.join(workingdir,file))
    



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
        number = numpy.random.randint(-1000000, -1)
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


    def test_return_eqaul_columns(self):
        """
        The assertion should areturn a number columns in retrieved dataset, should be equal to inserted dataset number of columns
        """

        # Random number of rows and columns
        r,c = numpy.random.randint(2,10,size=2)

        # Random array for dataset
        random_array = numpy.random.random([r,c])

        # Creating Dataset
        df = pandas.DataFrame(data = random_array,  columns = [str(i) for i in range(c)])

        # Inserting and retrieving
        result = self.insert_retrieve_dataframe('test_eqaul_columns','unittest_parquet',df)

        self.assertEqual(c,len(result.columns),'Number of rows not eqaul')

    
    def test_return_random_row_column_value(self):
        """
        The assertion should areturn a value from a random row and column in retrieved dataset, value should be equal to inserted dataset at same location
        """

        # Random number of rows and columns
        r,c = numpy.random.randint(2,10,size=2)

        # Random row and column to extract tested value
        r1 = numpy.random.randint(r)
        c1 = numpy.random.randint(c)

        # Random List for dataset
        random_array = numpy.random.random([r,c])

        # Creating Dataset
        df = pandas.DataFrame(data = random_array, columns = [str(i) for i in range(c)])
        
        # Sample value
        random_value = df.loc[r1,str(c1)]

        # Inserting and retrieving
        result = self.insert_retrieve_dataframe('test_random_value_in_df','unittest_parquet',df)

        self.assertEqual(random_value, result.loc[r1,str(c1)],'Random value not eqaul')
    
    def test_insert_value_overwrite_and_retrieve(self):
        """
        The assertion will be true if a dataset is stored in parquet and metadata,
        aftewords the dataset and metadata are overwritten and retrieved the expected updated value
        """
        overwrite_value = 11
        overwrite_index = 2
        df = pandas.DataFrame({'x':[0,1,2,3,4]})
        #inserting dataframe
        result = self.insert_retrieve_dataframe('test_overwrite_value_in_df','unittest_parquet',df)
        result.loc[overwrite_index,'x'] = overwrite_value
        result_after_overwrite = self.overwrite_parquet('test_overwrite_value_in_df','unittest_parquet',result)
        
        self.assertEqual(overwrite_value,result_after_overwrite.loc[overwrite_index,'x'],'Value after overwrite not equal')


    def test_update_metadata_and_retrieve(self):
        """
        The test modfies the metadata, afterwords the test reads the path to the parquet, the test will be successfull
        if the read path returns same value after metadata was modified and queried correctly based on path value to parquet file
        """
        number = numpy.random.randint(0,10000)
        self.insert_retrieve('test_update_metadata_1','unittest_parquet',number)
        item_to_update = datalayer.Measurements.getDocuments(projectName='unittest_parquet',test_id='test_update_metadata_1')[0]
        item_to_update.desc['test_id'] = "test_update_metadata_2"
        item_to_update.save()
        
        result = datalayer.Measurements.getDocuments(projectName='unittest_parquet',test_id='test_update_metadata_2')
        dataset_read = pandas.read_parquet(result[0].resource)
        result_val = dataset_read.loc[0,'x']
        self.assertEqual(number,result_val,'Return number not equal to given number after metadata update')

if __name__ =='__main__':
    unittest.main(warnings='ignore')
