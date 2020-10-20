import scipy.optimize as optimization 
import matplotlib.pyplot as plt 
import numpy
from scipy.stats import lognorm


class EstimateTBProbit(object): 
	"""
		This class helps to gusss the new sigma, as it is (for some odd reason) 
		does not converge. 
	"""
	def __init__(self,n):
		self.n = n 

	def transform(self,dosage): 
		return 2*(dosage/2)**self.n

	def plotGuess(self,newProbit,E50,probit): 
		"""
			Plots the transformed data (which is the correct) vs. 
			the guess of the lognormal parameter. 

			We recall that the transformation (assuming the exposure data is of 2 minutes) is: 
			
				TB-thresholds = 2*(Dosage/2)^n 

			Remember that everything is in base 10. 
		"""
		DoseList = numpy.arange(1,100,0.01) # all our dosages are < 100 . 
		gStd = 1/probit
		newgStd = 1/newProbit
		
		newE50 = self.transform(E50) 
		newDoseList = self.transform(DoseList) 
		a = numpy.log10(numpy.e)

		EffectList = lognorm.cdf(DoseList/a,s=gStd/a,loc=0,scale=E50/a)
		
		plt.semilogx(newDoseList,EffectList,label='The tenberge dose response') 

		guessEffectList = lognorm.cdf(newDoseList/a,s=newgStd/a,loc=0,scale=newE50/a) 
		plt.semilogx(newDoseList,guessEffectList,label='Guess')
		plt.legend()
		plt.show()
		print("New LD50 %s" % newE50)
		print("New gStd %s probit %s" % (probit/self.n,self.n/probit))

GBestimator = EstimateTBProbit(1.5)

##
# Civilians 
##

### Death 
# GB, LD50 18, probit = 4
#GBestimator.plotGuess(2.6666,18,4)

### Severe 
# GB, E50 13, probit = 4
#GBestimator.plotGuess(2.6666,13,4)

### Mild 
# GB, E50 8, probit = 3
#GBestimator.plotGuess(2.,8,3)

### Light
# GB, E50 0.5, probit = 3
#GBestimator.plotGuess(2,0.5,3)


##
# Soldiers
##

### Death 
# GB, LD50 35, probit = 12
#GBestimator.plotGuess(8,35,12)

### Severe 
# GB, E50 25, probit = 12
#GBestimator.plotGuess(8,25,12)

### Mild 
# GB, E50 16, probit = 11
#GBestimator.plotGuess(7.333,16,11)

### Light
# GB, E50 1, probit = 10
#GBestimator.plotGuess(6.666,1,10)


