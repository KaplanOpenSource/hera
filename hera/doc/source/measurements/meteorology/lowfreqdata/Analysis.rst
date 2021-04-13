Adding date and time details to the data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes we need to perform simple queries on the date. The
addDatesColumns function in the lowfreqdata.analysis module adds date
and time columns to the dataframe from the date index column. These
fields can, for example, when you want to extract data by different time
slices (by season, year, hour, etc.)

First, we get the data

.. code:: ipython3

    stationName = 'BET DAGAN'
    from hera.measurements.meteorology import lowfreqdata 
    datadb = lowfreqdata.lowfreqDataLayer.getStationDataFromDB(StationName=stationName)


.. parsed-literal::

    WARNING:root:FreeCAD not Found, cannot convert to STL


.. parsed-literal::

    You must install python-wrf to use this package 


.. code:: ipython3

    datadb.getData().compute().head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>stn_name</th>
          <th>stn_num</th>
          <th>BP</th>
          <th>DiffR</th>
          <th>Grad</th>
          <th>NIP</th>
          <th>RH</th>
          <th>Rain</th>
          <th>STDwd</th>
          <th>TD</th>
          <th>TDmax</th>
          <th>TDmin</th>
          <th>TG</th>
          <th>Time</th>
          <th>WD</th>
          <th>WDmax</th>
          <th>WS</th>
          <th>WS1mm</th>
          <th>WSmax</th>
          <th>Ws10mm</th>
        </tr>
        <tr>
          <th>time_obs</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>2020-09-01 00:00:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.4</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>7.9</td>
          <td>27.6</td>
          <td>27.6</td>
          <td>27.6</td>
          <td>27.4</td>
          <td>2351</td>
          <td>212</td>
          <td>223</td>
          <td>2.0</td>
          <td>2.2</td>
          <td>2.7</td>
          <td>2.0</td>
        </tr>
        <tr>
          <th>2020-09-01 00:10:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.3</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>8.7</td>
          <td>27.7</td>
          <td>27.8</td>
          <td>27.6</td>
          <td>27.9</td>
          <td>1</td>
          <td>208</td>
          <td>195</td>
          <td>1.9</td>
          <td>2.8</td>
          <td>3.3</td>
          <td>2.0</td>
        </tr>
        <tr>
          <th>2020-09-01 00:20:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.2</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>7.8</td>
          <td>27.9</td>
          <td>27.9</td>
          <td>27.8</td>
          <td>27.7</td>
          <td>18</td>
          <td>200</td>
          <td>191</td>
          <td>2.1</td>
          <td>2.6</td>
          <td>3.0</td>
          <td>2.1</td>
        </tr>
        <tr>
          <th>2020-09-01 00:30:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.2</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>10.4</td>
          <td>27.9</td>
          <td>27.9</td>
          <td>27.8</td>
          <td>27.6</td>
          <td>21</td>
          <td>199</td>
          <td>192</td>
          <td>1.9</td>
          <td>2.6</td>
          <td>3.0</td>
          <td>2.1</td>
        </tr>
        <tr>
          <th>2020-09-01 00:40:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.2</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>9.0</td>
          <td>27.8</td>
          <td>27.8</td>
          <td>27.7</td>
          <td>27.4</td>
          <td>37</td>
          <td>202</td>
          <td>187</td>
          <td>1.8</td>
          <td>2.5</td>
          <td>3.0</td>
          <td>1.9</td>
        </tr>
      </tbody>
    </table>
    </div>



Now we import the analysis function

.. code:: ipython3

    lowfreqdata.analysis.addDatesColumns(datadb.getData()).compute().head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>stn_name</th>
          <th>stn_num</th>
          <th>BP</th>
          <th>DiffR</th>
          <th>Grad</th>
          <th>NIP</th>
          <th>RH</th>
          <th>Rain</th>
          <th>STDwd</th>
          <th>TD</th>
          <th>...</th>
          <th>WS</th>
          <th>WS1mm</th>
          <th>WSmax</th>
          <th>Ws10mm</th>
          <th>curdate</th>
          <th>yearonly</th>
          <th>monthonly</th>
          <th>dayonly</th>
          <th>timeonly</th>
          <th>season</th>
        </tr>
        <tr>
          <th>time_obs</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>2020-09-01 00:00:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.4</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>7.9</td>
          <td>27.6</td>
          <td>...</td>
          <td>2.0</td>
          <td>2.2</td>
          <td>2.7</td>
          <td>2.0</td>
          <td>2020-09-01 00:00:00</td>
          <td>2020</td>
          <td>9</td>
          <td>1</td>
          <td>00:00:00</td>
          <td>Autumn</td>
        </tr>
        <tr>
          <th>2020-09-01 00:10:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.3</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>8.7</td>
          <td>27.7</td>
          <td>...</td>
          <td>1.9</td>
          <td>2.8</td>
          <td>3.3</td>
          <td>2.0</td>
          <td>2020-09-01 00:10:00</td>
          <td>2020</td>
          <td>9</td>
          <td>1</td>
          <td>00:10:00</td>
          <td>Autumn</td>
        </tr>
        <tr>
          <th>2020-09-01 00:20:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.2</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>7.8</td>
          <td>27.9</td>
          <td>...</td>
          <td>2.1</td>
          <td>2.6</td>
          <td>3.0</td>
          <td>2.1</td>
          <td>2020-09-01 00:20:00</td>
          <td>2020</td>
          <td>9</td>
          <td>1</td>
          <td>00:20:00</td>
          <td>Autumn</td>
        </tr>
        <tr>
          <th>2020-09-01 00:30:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.2</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>10.4</td>
          <td>27.9</td>
          <td>...</td>
          <td>1.9</td>
          <td>2.6</td>
          <td>3.0</td>
          <td>2.1</td>
          <td>2020-09-01 00:30:00</td>
          <td>2020</td>
          <td>9</td>
          <td>1</td>
          <td>00:30:00</td>
          <td>Autumn</td>
        </tr>
        <tr>
          <th>2020-09-01 00:40:00</th>
          <td>BET DAGAN</td>
          <td>54</td>
          <td>1001.2</td>
          <td>0</td>
          <td>-9999</td>
          <td>-9999</td>
          <td>0</td>
          <td>0</td>
          <td>9.0</td>
          <td>27.8</td>
          <td>...</td>
          <td>1.8</td>
          <td>2.5</td>
          <td>3.0</td>
          <td>1.9</td>
          <td>2020-09-01 00:40:00</td>
          <td>2020</td>
          <td>9</td>
          <td>1</td>
          <td>00:40:00</td>
          <td>Autumn</td>
        </tr>
      </tbody>
    </table>
    <p>5 rows × 26 columns</p>
    </div>



By default, the function takes the ‘index’ column as the source for
calendar information. The user can give his own date and month columns
using the ‘datecolumn’ and ‘monthcolumn’ variables. In addition, the
default function handles dask-type data. You can change the type to
‘pandas’ with the ‘dataType’ variable.

Calculating hourly distribution of a field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thid procedure calculates the hourly distribution an field for a
specified number of bin.

.. code:: ipython3

    x_mid,y_mid,M = lowfreqdata.analysis.calcHourlyDist(datadb.getData(),Field='WS')

