class errorCalculation():

    def calculateFB(self,data, modelColumn='model', measureColumn='measure'):
        """
        Calculates the fractional mean bias.

        :param data: pandas of the model and raw data together.
        :param modelColumn: The name of the model's data column.
        :param measureColumn: The name of the experimental measured data.
        :return: float FB.
        """

        FB = 2*(data[modelColumn]-data[measureColumn]).mean()/(data[modelColumn].mean()+data[measureColumn].mean())
        return FB


    def calculateNMSE(self,data, modelColumn='model', measureColumn='measure'):
        """
        Calculates the normalized mean-square error.

        :param data: pandas of the model and raw data together.
        :param modelColumn: The name of the model's data column.
        :param measureColumn: The name of the experimental measured data.
        :return: float NMSE.
        """

        NMSE = ((data[modelColumn]-data[measureColumn])**2).mean()/(data[modelColumn].mean()*data[measureColumn].mean())
        return NMSE


    def calculateFAC(self,data, relation=2, modelColumn='model', measureColumn='measure'):
        """
        Calculates the FAC criteria.

        :param data: pandas of the model and raw data together.
        :param modelColumn: The name of the model's data column.
        :param measureColumn: The name of the experimental measured data.
        :return: float FAC2.
        """

        model_over_measure = data[modelColumn]/data[measureColumn]
        FAC2 = model_over_measure.apply(lambda x: 1/relation<x<relation).sum()/model_over_measure.count()
        return FAC2


    def calculateNAD(self,data, modelColumn='model', measureColumn='measure'):
        """
        Calculates the normalized absolute difference.

        :param data: pandas of the model and raw data together.
        :param modelColumn: The name of the model's data column.
        :param measureColumn: The name of the experimental measured data.
        :return: float NAD.
        """

        NAD = abs((data[modelColumn]-data[measureColumn]).mean())/(data[modelColumn].mean()+data[measureColumn].mean())
        return NAD


    def calculateR(self,data, modelColumn='model', measureColumn='measure'):
        """
        Calculates the linear correlation.

        :param data: pandas of the model and raw data together.
        :param modelColumn: The name of the model's data column.
        :param measureColumn: The name of the experimental measured data.
        :return: float R.
        """

        R = ((data[modelColumn]-data[modelColumn].mean())*(data[measureColumn]-data[measureColumn].mean())).mean()/(data[modelColumn].std()*data[measureColumn].std())
        return R


    def calculateFB_FN(self,data, modelColumn='model', measureColumn='measure'):
        """
        Calculates FB_FN that tell us about under estimate in FB criteria.

        :param data: pandas of the model and raw data together.
        :param modelColumn: The name of the model's data column.
        :param measureColumn: The name of the experimental measured data.
        :return: float FB_FN.
        """

        FB_FN = sum(abs(data[modelColumn]-data[measureColumn])+(data[modelColumn]-data[measureColumn]))/sum(data[modelColumn]+data[measureColumn])
        return FB_FN


    def calculateFB_FP(self,data, modelColumn='model', measureColumn='measure'):
        """
        Calculates FB_FP that tell us about over estimate in FB criteria.

        :param data: pandas of the model and raw data together.
        :param modelColumn: The name of the model's data column.
        :param measureColumn: The name of the experimental measured data.
        :return: float FB_FP.
        """

        FB_FP = sum(abs(data[modelColumn]-data[measureColumn])+(data[measureColumn]-data[modelColumn]))/sum(data[modelColumn]+data[measureColumn])
        return FB_FP