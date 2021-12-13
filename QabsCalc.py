#!/usr/bin/python
import numpy

'''
Description of code here.
'''

def QabsCalc(grainsize, lambda, Qabs, lambdaQabs, a):

# The dust sizes are stored in a, wavelengths in lambdaQabs, & the absorption 
# coefficients in Qabs
# Qabs is only calculated for grain sizes larger than 0.01 microns & smaller 
# than 3.0 microns, if the user supplied grain size is smaller than 0.01, 
# the program uses the Qabs for 0.01.  If the grain size is larger than 3 
# microns, it is approximated as a power law.  This subroutine searches a list
# of precalculated absorption coefficients for 21 different grain sizes & 112
# different wavelengths.  If the specified grain size & wavelength are not on
# the list then it interpolates between the next lowest & highest value.
   lambda0=lambda*1e-4
   index=-1

   if (grainsize >= 0.01) AND (grainsize <= 3): 
	  index=(a EQ grainsize).nonzero()	#
	  index = index[0]
   if (ARRAY_EQUAL(index,-1)):
		uindex = (a > grainsize).nonzero()	#FIX
		uindex = uindex[0]
		upperindex=uindex[0]
		lindex=(a < grainsize).nonzero()	#FIX
		lindex = lindex[0]
		lowerindex=size(lindex)
		lowerindex=lowerindex(1)-1 

   if ((lambda >= 0.01) AND (lambda <= 36000)):
	   lambdaindex=(lambda == lambdaQabs).nonzero()	#FIX
	   lambdaindex = lambdaindex[0]
	  if (ARRAY_EQUAL(lambdaindex,-1)):
		 uLambdaIndex=(lambdaQabs > lambda).nonzero()	#FIX
		 uLambdaIndex = uLambdaIndex[0]
		 upperLambdaIndex=uLambdaIndex[0]
		 LambdaIndex=(lambdaQabs < lambda).nonzero()	#FIX
		 LambdaIndex = lLambdaIndex[0]
		 lowerLambdaIndex=len(lLambdaIndex)
		 lowerLambdaIndex=lowerLambdaIndex(1)-1
		 if (index[0] == -1):
			f=[Qabs(lowerindex,lowerLambdaIndex),Qabs(upperindex,upperLambdaIndex)]
			Qabsreturn=interpolate(f,lambda)
		 else:
			f=[Qabs(index[0],lowerLambdaIndex),Qabs(index[0],upperLambdaIndex)]
			Qabsreturn=interpolate(f,lambda)
	  else:   
		 if (index[0] EQ -1):
			f=[Qabs(lowerindex,lambdaindex),Qabs(upperindex,lambdaindex)]
			Qabsreturn=interpolate(f,lambda)
		 else:
			Qabsreturn=Qabs(index[0],lambdaindex)
   else:
	  if (lambda <= 0.01):
		 lambdaindex=0
		 if (lambda GE 36000): 
			lambdaindex=111

   if ((lambda <= 0.01) OR (lambda >= 36000)):
	  if (index[0] == -1): 
		 if (grainsize < 0.01): 
			Qabsreturn=Qabs(0,lambdaindex)
		 else:
			f=[Qabs(lowerindex,lambdaindex),Qabs(upperindex,lambdaindex)]
			Qabsreturn=interpolate(f,lambda)
	  else:
		 Qabsreturn=Qabs(index[0],lambdaindex)
return Qabsreturn