#!/usr/bin/python
import scipy

'''
# The dust sizes are stored in a, wavelengths in lambdaQabs, & the absorption 
# coefficients in Qabs
# Qabs is only calculated for grain sizes larger than 0.01 microns & smaller 
# than 3.0 microns, if the user supplied grain size is smaller than 0.01, 
# the program uses the Qabs for 0.01.  If the grain size is larger than 3 
# microns, it is approximated as a power law.  This subroutine searches a list
# of precalculated absorption coefficients for 21 different grain sizes & 112
# different wavelengths.  If the specified grain size & wavelength are not on
# the list then it interpolates between the next lowest & highest value.

CORBIN: The main problems are with scipy.interpolate() and using operators on booleans 
'''

def QabsCalc(grainsize, lambda_in, Qabs, lambdaQabs, a):
	lambda0=lambda_in*1e-4
	index=-1

	if grainsize >= 0.01 and grainsize <= 3:
		a = grainsize
		index_o: object=(a)
		index = index_o[0]
	if index==-1:
		uindex = (a > grainsize)
		uindex = uindex[0]
		upperindex=uindex[0]
		lindex=(a < grainsize)
		lindex = lindex[0]
		lowerindex=len(lindex)
		lowerindex=lowerindex-1

	if lambda_in >= 0.01 and lambda_in <= 36000:
		lambdaindex=(lambda_in == lambdaQabs)
		lambdaindex = lambdaindex[0]
		if lambdaindex==-1:
			uLambdaIndex=(lambdaQabs > lambda_in)
			uLambdaIndex = uLambdaIndex[0]
			upperLambdaIndex=uLambdaIndex[0]
			lLambdaIndex=(lambdaQabs < lambda_in)
			lLambdaIndex = lLambdaIndex[0]
			lowerLambdaIndex=len(lLambdaIndex)
			lowerLambdaIndex=lowerLambdaIndex-1
			if index[0] == -1:
				f=[Qabs(lowerindex,lowerLambdaIndex),Qabs(upperindex,upperLambdaIndex)]
				Qabsreturn=scipy.interpolate(f,lambda_in)
			else:
				f=[Qabs(index[0],lowerLambdaIndex),Qabs(index[0],upperLambdaIndex)]
				Qabsreturn=scipy.interpolate(f,lambda_in)
		else:
			if (index[0]==-1):
				f=[Qabs(lowerindex,lambdaindex),Qabs(upperindex,lambdaindex)]
				Qabsreturn=scipy.interpolate(f,lambda_in)
			else:
				Qabsreturn=Qabs(index[0],lambdaindex)
	else:
		if lambda_in <= 0.01:
			lambdaindex=0
		if lambda_in >= 36000:
			lambdaindex=111

	if ((lambda_in <= 0.01 or lambda_in >= 36000)):
		if (index[0] == -1):
			if (grainsize < 0.01):
				Qabsreturn=Qabs(0,lambdaindex)
			else:
				f=[Qabs(lowerindex,lambdaindex),Qabs(upperindex,lambdaindex)]
				Qabsreturn=scipy.interpolate(f,lambda_in)
	else:
		Qabsreturn=Qabs(index[0],lambdaindex)

	return Qabsreturn