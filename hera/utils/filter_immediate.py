class Filter:
    """Apply threshold and interval filters to a pandas DataFrame."""

    def __init__(self, data, inplace = False):
        """
        Initialize the Filter.

        Parameters
        ----------
        data : DataFrame
            The data to filter.
        inplace : bool
            If True, mutate this instance; otherwise return a new Filter.
        """
        self.inplace = inplace
        self.data = data

    def threshold(self, preposition, bound, column = None):
        """
        Filter rows by a comparison operator on a column or the index.

        Parameters
        ----------
        preposition : str
            One of lt, lte, gt, gte, abs_lt, abs_lte, abs_gt, abs_gte, eq, neq.
        bound : scalar
            The threshold value.
        column : str, optional
            Column to filter on; uses the index if None.
        """
        data_column = (self.data.index if column is None else self.data[column])

        if preposition == "lt":
            data_masked = self.data.loc[data_column < bound]
        elif preposition == "lte":
            data_masked = self.data.loc[data_column <= bound]
        elif preposition == "gt":
            data_masked = self.data.loc[data_column > bound]
        elif preposition == "gte":
            data_masked = self.data.loc[data_column >= bound]
        elif preposition == "abs_lt":
            data_masked = self.data.loc[abs(data_column) < bound]
        elif preposition == "abs_lte":
            data_masked = self.data.loc[abs(data_column) <= bound]
        elif preposition == "abs_gt":
            data_masked = self.data.loc[abs(data_column) > bound]
        elif preposition == "abs_gte":
            data_masked = self.data.loc[abs(data_column) >= bound]
        elif preposition == "eq":
            data_masked = self.data.loc[data_column == bound]
        elif preposition == "neq":
            data_masked = self.data.loc[data_column != bound]
        else:
            raise ValueError("preposition must be one of the recognized prepositions: lt,lte,gt,gte,abs_lt,abs_lte,abs_gt,abs_gte,eq,neq")

        if self.inplace:
            self.data = data_masked
            return self
        else:
            return Filter(data = data_masked, inplace = self.inplace)


    def outsideInterval(self,  lower_bound, upper_bound, column = None):
        """
        Keep rows whose values fall outside [lower_bound, upper_bound).

        Parameters
        ----------
        lower_bound : scalar
            Lower bound of the excluded interval.
        upper_bound : scalar
            Upper bound of the excluded interval.
        column : str, optional
            Column to filter on; uses the index if None.
        """
        data_column = (self.data.index if column is None else self.data[column])

        data_masked = self.data.loc[(data_column < lower_bound) | (data_column >= upper_bound)]

        if self.inplace:
            self.data = data_masked
            return self
        else:
            return Filter(data = data_masked, inplace = self.inplace)