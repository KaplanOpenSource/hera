import pandas
from argos.manager import  DEPLOY, DESIGN
from hera.datalayer import datatypes

class experimentAnalysis:
    """
        A basic class for the analysis of the data.
    """
    TECHNICALDOC_FREQUENCY = "technical_doc_frequency"  # A document type for caching of the frequency
    _datalayer = None


    @property
    def datalayer(self):
        return self._datalayer

    def __init__(self,datalayer):
        self._datalayer = datalayer

    def getDeviceLocations(self,entityTypeName,trialName,trialSetName=None,trialState=DEPLOY):
        """
        Returns a pandas with the device locations.

        Parameters
        ----------
        entityTypeName: str
                The entity type.

        trialName: str
                The trial name.

        trialSetName: str
                The trial set name.

        trialState: str
            DEPLOY or DESIGN.

        Returns
        -------
            pandas.DataFrame
        """
        trialSetName = self.datalayer.defaultTrialSet if trialSetName is None else trialSetName
        if trialSetName is None:
            raise ValueError("trialSetName is None. Either set default with the getExperiment, or set explicitly here")
        return self.datalayer.trialSet[trialSetName][trialName].entitiesTable(trialState).query("entityType==@entityTypeName")


    def getTurbulenceStatistics(self,sonicData,samplingWindow,height=1):
        """
        Returns the turbulence analysis of the sonic/kaijo data for the start->end times.

        Parameters
        ----------
        sonicData: pandas/Dask
            the data of the sonic to analyze.

        samplingWindow: str
            The window width for analysis

        height: int
            Height

        Returns
        -------
            singlePointTurbulenceStatistics class
        """

        highfreqtk = self.datalayer.toolkitExtension.sonicHighFreqToolkit
        analysis = highfreqtk.analysis.singlePointTurbulenceStatistics(sonicData=kaijoData,
                                                                       start=None,
                                                                       end=None,
                                                                       height=kaijoHeight,
                                                                       samplingWindow=samplingWindow,
                                                                       buildingHeight= 0,
                                                                       averagedHeight=0,
                                                                       isMissingData=True)

        return analysis

    def _splitName(self,x):
        if len(x.split(" ")) > 1:
            return x.split(" ")[1]
        else:
            return x.split(" ")[0]

    def getDeviceTypeTransmissionFrequencyOfTrial(self,
                                                  deviceType,
                                                  trialName,
                                                  trialSetName=None,
                                                  samplingWindow="1min",
                                                  normalize=True,
                                                  completeTimeSeries=True,
                                                  completeDevices=True,
                                                  wideFormat=True,
                                                  recalculate=True):
        """
        Returns the frequency of a device in the requested trialSet/trial in a wide format (default) or long format (wideFormat=False).
        The normalized frequency is obtained from the table defined in this class.
        Normalize to the optimal sampling rate if normalize is true.
        Can also complete the devices that were planned but did not transmit at all (and therefore, do not appear in the DB)
        or the time within the trial were no transmission was obtained.

        Parameters
        ----------
        deviceType : string
                The type of the device to present.

        trialName : string
                The name of the trial to show.

        trialSetName : string
                The name of the trial set.

        samplingWindow : string
                A time string (of pandas). i.e. '1min' and ect.

        normalize : bool
                If true, normalize to the planned message rate

        completeTimeSeries : bool
                If true, and the release has end and start times, then add all the missing time slots
                with the requested interval.
                default: True

        completeDevices : bool=True
                If True, add all the missing devices that were related to the trial.
                default: True

        wideFormat: bool
                If true, return in wide format, else return in long format
                default = wide format

        recalculate: bool
                Recalculate the cache of the data.
                default False.

        Returns
        ----------
            pandas.DataFrame in wide or long formats.
        """

        trialSetName = self.datalayer.defaultTrialSet if trialSetName is None else trialSetName
        experimentData = self.datalayer

        qry = dict(deviceType = deviceType,
                   samplingWindow = samplingWindow,
                   trialName = trialName,
                   trialSetName = trialSetName,
                   completeTimeSeries=completeTimeSeries,
                   completeDevices=completeDevices)

        docList = experimentData.getCacheDocuments(type=self.TECHNICALDOC_FREQUENCY,**qry)

        if len(docList) ==0 or recalculate:
            # calculating the cache and updating it / or writing a new one.
            data = experimentData.trialSet[trialSetName][trialName].getData(deviceType)
            colCount = data.columns[1]
            freq = data.groupby('deviceName').resample("1min").count()[colCount].reset_index()

            freq = freq.assign(deviceID=freq.deviceName.apply(lambda x: self._splitName(x))).sort_values("deviceID").rename(columns=({colCount:"Frequency"}))

            freq.set_index('timestamp')
            pvt = freq.pivot(index="timestamp", columns="deviceName", values="Frequency")
            if completeDevices:
                for deviceName in experimentData.entityType[deviceType].keys():
                    if deviceName not in pvt:
                        pvt[deviceName] = 0

            # Filling in time (only if trial has start and end).
            tstart = experimentData.trialSet[trialSetName].trials.query("trialName == @trialName ")['TrialStart'].iloc[0]
            tend = experimentData.trialSet[trialSetName].trials.query("trialName == @trialName ")['TrialEnd'].iloc[0]

            if (tend is not None and completeTimeSeries):
                actualDates = pandas.date_range(tstart, tend, freq="1min")

                missingTimes = actualDates[[x not in pvt.index for x in actualDates]]
                for tme in missingTimes:
                    pvt.loc[tme, :] = 0

            pvt = pvt.sort_index(ascending=False)

            # update the DB.
            if len(docList) == 0 :
                # new record
                experimentData.addCacheDocument(type=self.TECHNICALDOC_FREQUENCY,
                                                  dataFormat=datatypes.JSON_PANDAS,
                                                  resource=pvt.to_json(),
                                                  desc=qry)

            else:
                # update existing one
                doc = docList[0]
                doc.resource = pvt.to_json()
                doc.save()
        else:
            # return the data from the cache.
            pvt = docList[0].getData()

        if normalize:
            total_expected_messages = self.getDeviceTypePlannedMessageCount(deviceType=deviceType,
                                                                            samplingWindow=samplingWindow)
            pvt = pvt / total_expected_messages


        if not wideFormat:
            ret = pandas.melt(frame=pvt.reset_index(),id_vars="index",value_vars=pvt.columns,var_name="deviceName",value_name="Frequency").rename(columns=dict(index="timestamp"))
        else:
            ret = pvt

        return ret


    def getDeviceTypePlannedMessageCount(self,deviceType,samplingWindow="1min"):
        """
        Returns the number of messages that should be obtained in a sampling window,according to the device type.

        Parameters
        ----------
        deviceType : string
                The type of the device to present.

        samplingWindow : string
                A time string (of pandas). i.e. '1min' and ect.

        """



        samplingWindow_seconds = pandas.to_timedelta(samplingWindow).total_seconds()
        total_expected_messages = self.getOptimalFrequencyHz(deviceType)*samplingWindow_seconds

        if total_expected_messages<1:
            raise ValueError(f"Sampling window is too short. Use sampling window of at least {1/self.getOptimalFrequencyHz(deviceType)}s")

        return total_expected_messages





    def addMetadata(self,dataset,trialName,trialState=DEPLOY,trialSetName=None):
        """
        Adds the trial metadata to the dataset by joining on the device name.

        Parameters
        ----------
        dataset :  pandas/dask
                The dataset to merge. Assumes there is a field deviceName
        trialName : str
                The name of the trial to join

        trialState : str
                The state of the trial (Design/Deploy).

        trialSetName : str [optional, default None]
                The name of hte trialset. If None, use the default trialset.

        Returns
        -------
            pandas/dask with the metadata.
        """
        trialSetName =  self.datalayer.defaultTrialSet if  trialSetName is None else trialSetName

        devicemetadata = self.datalayer.trialSet[trialSetName][trialName].entitiesTable(trialState)
        if len(devicemetadata) > 0:
            ret = dataset.reset_index().merge(devicemetadata, left_on="deviceName", right_on="entityName").set_index("timestamp")
        else:
            ret=None

        return ret

    def addTrialProperties(self,data,trialName,trialSetName=None):
        """
        Adds the calculated properties of the trial to the data.
        - FromRelease: seconds from the release.
        - FromStart  : seconds from the start.
        Assume the time column is timestamp.

        Parameters
        ----------
        data: pandas/dask
            The data.
        trialName: str
            The trial name.
        trialSetName:str
            The trial set name.

        Returns
        -------
            pandas/dask
        """

        trialSetName = self.datalayer.defaultTrialSet if trialSetName is None else trialSetName

        startTime = self.datalayer.trialSet[trialSetName][trialName].properties['TrialStart']
        releaseTime = self.datalayer.trialSet[trialSetName][trialName].properties.get("ReleaseStart",None)

        if 'timestamp' not in data:
            tmp = data.reset_index()
            if 'timestamp' not in tmp:
                raise ValueError("The column 'timestamp' is not part of data (or in its index)")
        else:
            tmp = data

        tmp = tmp.assign(fromStart = tmp.timestamp-startTime,
                         fromRelease=tmp.timestamp-releaseTime)

        tmp = tmp.assign(fromStartSeconds = tmp.fromStart.apply(lambda x: x.total_seconds()),
                         fromReleaseSeconds=tmp.fromRelease.apply(lambda x: x.total_seconds()))


        return tmp




