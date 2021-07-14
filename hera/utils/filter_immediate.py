class Filter:
    def __init__(self, data, inplace = False):
        self.inplace = inplace
        self.data = data

    def threshold(self, preposition, bound, column = None):

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

        data_column = (self.data.index if column is None else self.data[column])

        data_masked = self.data.loc[(data_column < lower_bound) | (data_column >= upper_bound)]

        if self.inplace:
            self.data = data_masked
            return self
        else:
            return Filter(data = data_masked, inplace = self.inplace)