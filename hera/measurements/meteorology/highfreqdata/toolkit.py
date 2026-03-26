import os
import warnings

from .analysis.analysislayer import RawdataAnalysis
from hera import toolkit
from hera.datalayer import datatypes
from hera.utils.logging import get_classMethod_logger
from .parsers.CampbellBinary import Parser
from .parsers.TOA5 import ASCIIParser
from .parsers import normalise_sonic_df, normalise_trh_df, detect_device_type


class HighFreqToolKit(toolkit.abstractToolkit):
    """Toolkit for high-frequency (10–20 Hz) sonic anemometer and TRH data.

    Manages parsing raw binary/ASCII sensor files, normalising the output,
    storing data as versioned data sources in the project, and providing
    an analysis layer for turbulence statistics.

    Supported raw formats:

    - **Campbell Scientific TOB1 binary** — via :class:`Parser`
    - **Campbell Scientific TOA5 ASCII** — via :class:`ASCIIParser`
    """

    DOCTYPE_STATIONS = 'StationsData'
    DOCTYPE_MEASUREMENTS = 'MeasurementsData'

    def __init__(self, projectName, filesDirectory=None, connectionName=None):
        """Initialise the high-frequency meteorology toolkit.

        Parameters
        ----------
        projectName : str
            The project name.
        filesDirectory : str, optional
            Directory for saving parquet files. If ``None``, uses the
            project default.
        connectionName : str, optional
            Database connection name.
        """
        super().__init__(
            projectName=projectName,
            toolkitName="highFreqMeteorology",
            filesDirectory=filesDirectory,
            connectionName=connectionName,
        )
        logger = get_classMethod_logger(self, "init")
        logger.info("Init High frequency data")
        self._analysis = RawdataAnalysis(self)

    @property
    def docType(self):
        """str : Document type identifier for high-frequency data."""
        return f"{self.toolkitName}_HighFreqData"

    # ------------------------------------------------------------------
    # Parsing + normalising
    # ------------------------------------------------------------------

    def parseData(self, path, fromTime=None, toTime=None, parser="auto"):
        """Parse a raw data file and return normalised DataFrames.

        Automatically detects the file format (Campbell binary vs TOA5
        ASCII) unless *parser* is specified explicitly. The output is
        normalised to the standard column format expected by the analysis
        layer (lowercase ``u, v, w, T`` with a ``DatetimeIndex``).

        Parameters
        ----------
        path : str
            Path to a raw data file or a directory of files.
        fromTime : str or None, optional
            Start time filter.
        toTime : str or None, optional
            End time filter.
        parser : str, optional
            Parser to use: ``'auto'`` (default), ``'campbell'``, or
            ``'toa5'``.

        Returns
        -------
        list of tuple(pandas.DataFrame, dict)
            Each element is ``(normalised_dataframe, metadata_dict)``.
            A single-device file produces a list with one element.
            Multi-device files (TOA5) produce one element per device.
        """
        raw = self._parse_raw(path, fromTime, toTime, parser)
        return self._normalise_raw(raw)

    def loadData(self, name, path, fromTime=None, toTime=None,
                 parser="auto", version=(0, 0, 1), overwrite=False,
                 metadata=None):
        """Parse a raw data file, save as parquet, and register as a data source.

        This is the recommended way to ingest new high-frequency data. It
        parses, normalises, writes a parquet file to the toolkit's file
        directory, and registers the result as a versioned data source in
        the project.

        Parameters
        ----------
        name : str
            Data source name (e.g. ``'sonic_10m'``). For multi-device
            files, device names are appended (``'sonic_10m_Raw_Sonic_1'``).
        path : str
            Path to a raw data file or directory.
        fromTime : str or None, optional
            Start time filter.
        toTime : str or None, optional
            End time filter.
        parser : str, optional
            ``'auto'``, ``'campbell'``, or ``'toa5'``.
        version : tuple of int, optional
            Data source version ``(major, minor, patch)``. Default ``(0, 0, 1)``.
        overwrite : bool, optional
            If ``True``, overwrite an existing data source with the same name.
        metadata : dict, optional
            Additional metadata to store with the data source.

        Returns
        -------
        list of documents
            The created data source document(s).
        """
        results = self.parseData(path, fromTime, toTime, parser)

        # Ensure output directory exists
        parquet_dir = os.path.join(self.filesDirectory, "parquet")
        os.makedirs(parquet_dir, exist_ok=True)

        docs = []
        for df, meta in results:
            # Build data source name: append device info for multi-device
            if len(results) == 1:
                ds_name = name
            else:
                suffix = meta.get("deviceName", meta.get("deviceType", ""))
                ds_name = f"{name}_{suffix}" if suffix else name

            # Merge user metadata
            desc = dict(meta)
            if metadata:
                desc.update(metadata)

            # Save parquet
            output_file = os.path.join(parquet_dir, f"{ds_name}.parquet")
            df.to_parquet(output_file)

            # Register as data source
            doc = self.addDataSource(
                dataSourceName=ds_name,
                resource=output_file,
                dataFormat=datatypes.PARQUET,
                version=list(version),
                overwrite=overwrite,
                **desc,
            )
            docs.append(doc)

        return docs if len(docs) > 1 else docs[0]

    # ------------------------------------------------------------------
    # Legacy methods (backward compatible, now normalise output)
    # ------------------------------------------------------------------

    def campbelToParquet(self, binaryFile, fromTime=None, toTime=None,
                         chunkSize=10000):
        """Parse a Campbell binary file and return a normalised DataFrame.

        .. deprecated::
            Use :meth:`parseData` or :meth:`loadData` instead.

        Parameters
        ----------
        binaryFile : str
            Path to the Campbell binary ``.dat`` file.
        fromTime : str or None, optional
            Start time filter.
        toTime : str or None, optional
            End time filter.
        chunkSize : int, optional
            Records per batch (default 10000).

        Returns
        -------
        pandas.DataFrame
            Normalised DataFrame with lowercase columns and DatetimeIndex.
        """
        warnings.warn(
            "campbelToParquet() is deprecated. Use parseData() or loadData() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        campelParser = Parser(chunkSize=chunkSize)
        raw_df = campelParser.parse(path=binaryFile, fromTime=fromTime, toTime=toTime)
        results = self._normalise_raw(raw_df)
        if len(results) == 1:
            return results[0][0]
        return {meta.get("deviceName", f"device_{i}"): df for i, (df, meta) in enumerate(results)}

    def asciiToParquet(self, path, fromTime=None, toTime=None):
        """Parse a TOA5 ASCII file and return normalised DataFrames.

        .. deprecated::
            Use :meth:`parseData` or :meth:`loadData` instead.

        Parameters
        ----------
        path : str
            Path to a ``.dat`` file or directory.
        fromTime : str or None, optional
            Start time filter.
        toTime : str or None, optional
            End time filter.

        Returns
        -------
        dict
            Mapping of device names to normalised DataFrames.
        """
        warnings.warn(
            "asciiToParquet() is deprecated. Use parseData() or loadData() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        asciiParser = ASCIIParser()
        raw = asciiParser.parse(path=path, fromTime=fromTime, toTime=toTime)
        results = self._normalise_raw(raw)
        return {meta.get("deviceName", f"device_{i}"): df for i, (df, meta) in enumerate(results)}

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _parse_raw(self, path, fromTime, toTime, parser):
        """Parse raw data using the specified or auto-detected parser.

        Returns the raw (un-normalised) parser output.
        """
        if parser == "auto":
            parser = self._detect_parser(path)

        if parser == "campbell":
            p = Parser()
            return p.parse(path=path, fromTime=fromTime, toTime=toTime)
        elif parser == "toa5":
            p = ASCIIParser()
            return p.parse(path=path, fromTime=fromTime, toTime=toTime)
        else:
            raise ValueError(f"Unknown parser '{parser}'. Use 'auto', 'campbell', or 'toa5'.")

    @staticmethod
    def _detect_parser(path):
        """Auto-detect parser type from file content.

        Checks the first bytes of the file: TOB1 binary files start with
        a text header line containing ``'TOB1'``; TOA5 files start with
        ``'TOA5'`` or are CSV-like. Falls back to ``'campbell'`` for
        ``.dat`` files.
        """
        target = path
        if os.path.isdir(path):
            # Pick first .dat file in the directory
            for f in os.listdir(path):
                if f.endswith(".dat"):
                    target = os.path.join(path, f)
                    break
            else:
                raise FileNotFoundError(f"No .dat files found in {path}")

        try:
            with open(target, "rb") as fh:
                header = fh.read(512)
            # TOA5 ASCII files have 'TOA5' near the beginning
            if b"TOA5" in header[:100]:
                return "toa5"
            # TOB1 binary files have 'TOB1' near the beginning
            if b"TOB1" in header[:100]:
                return "campbell"
        except Exception:
            pass

        # Fallback: try to read as CSV (TOA5)
        try:
            with open(target, "r") as fh:
                first_line = fh.readline()
            if "," in first_line and any(kw in first_line for kw in ("Raw_Sonic", "Tct_Trh", "TIMESTAMP")):
                return "toa5"
        except Exception:
            pass

        # Default to campbell for binary .dat files
        return "campbell"

    @staticmethod
    def _normalise_raw(raw):
        """Normalise raw parser output to a list of (DataFrame, metadata).

        Handles both single-DataFrame (CampbellBinary) and dict-of-DataFrames
        (TOA5) output.
        """
        import pandas as pd

        if isinstance(raw, dict):
            # TOA5 returns dict[device_name → DataFrame]
            results = []
            for device_name, df in raw.items():
                dtype = detect_device_type(df)
                if dtype == "sonic":
                    norm_df, meta = normalise_sonic_df(df)
                elif dtype == "TRH":
                    norm_df, meta = normalise_trh_df(df)
                else:
                    # Unknown type — try sonic normalisation
                    norm_df, meta = normalise_sonic_df(df)
                meta["deviceName"] = device_name
                results.append((norm_df, meta))
            return results

        elif isinstance(raw, pd.DataFrame):
            # CampbellBinary returns a single DataFrame
            dtype = detect_device_type(raw)
            if dtype == "sonic":
                norm_df, meta = normalise_sonic_df(raw)
            elif dtype == "TRH":
                norm_df, meta = normalise_trh_df(raw)
            else:
                norm_df, meta = normalise_sonic_df(raw)
            return [(norm_df, meta)]

        else:
            raise TypeError(f"Unexpected parser output type: {type(raw)}")
