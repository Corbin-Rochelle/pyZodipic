#!/usr/bin/python
import numpy
from scipy.ndimage.interpolation import shift

'''
CORBIN: This seems to be working as intended 
'''

def temperaturecalc(lstar, tstar, lambda0microns, rau, t, lambdaQabs, emit):
   # program to calculate the equilibrium
   # temperature (t) of
   # dust grains of arbitrary size, in
   # microns (lambda0microns)
   # as a function of the distance from
   # the star in AU, (rau)
   
   # calculate stellar radius from effective temperature 
   # & luminosity
   # rstar**2*(tstar/5770.0)**4.0=lstar
    rstar=(lstar**0.5)*(tstar/5770.0)**(-2.0)  # solar units
   # radius of the star in cm
    rstarcm=rstar*6.9599e10
   # effective dust size, in microns
    lambda0=lambda0microns*10
   # tol is the fractional error in the final temperature
    tol=0.001
   #****************************************
   # fundamental constants
    kb=1.38066e-16
    c=2.9979e10
    h=6.62608e-27
   # ***************************************
   # Uses a stored list of 112 wavelengths to calculate nu.
    nu = c/(lambdaQabs*1e-4)  # cgs
    deltanu=abs(nu-shift(nu,-1))
    deltanu[111]=0.0
    nu[111]=0.0
    sz=len(rau)
    raulen=sz
    t=numpy.zeroes(raulen)
   
    for i in range(0, raulen):
   # distance from the star in cm
        rcm=rau[i]*1.495979e13
   # guess that the grain is at blackbody temperature 
        tgrain=278.0
   #  Calculate stellar Planck spectrum Bnu (erg s**-1 cm**-2 ster**-1 Hz**-1)
        xb=h*nu/(kb*tstar)
        bnu=xb**3.0/(numpy.exp(xb)-1.0)
        bnu=bnu*2.0*((kb*tstar)**3.0)/((c*h)**2.0)
        bnu[111]=0.0
	  
   # efficiency of absorption
   # If the dust size is greater than 3 microns, the absorption is approximated
   # as a power law, if it is less than 3 microns it uses the appropriate
   # absorption coefficient for each wavelength.
   
    eabsorb=((rstarcm/rcm)**2.0)*sum(bnu*deltanu*emit)
    ok=0
    bad=0
    tstep=20.0
    while (ok == 0):
   # calculate emitted energy
   #Calculate Planck spectrum Bnu (erg s**-1 cm**-2 ster**-1 Hz**-1)
        xb=h*nu/(kb*tgrain)
        bnu=xb**3.0/(numpy.exp(xb)-1.0)
        bnu=bnu*2.0*((kb*tgrain)**3.0)/((c*h)**2.0)
        bnu[111]=0.0
	   
   # If the dust size is greater than 3 microns, the emission is approximated
   # as a power law, if it is less than 3 microns it gets the appropriate
   # absorption coefficient for each wavelength from QabsCalc. 
        eemit=4.0*sum(bnu*deltanu*emit)
        match=eemit/eabsorb
        oldbad=bad
   
   # If the guess is too low, increase the temperature
        if match > 1.0+tol:
            bad=1
            tgrain=tgrain-tstep
   
   
   # If the guess is too high, decrease the temperature
        if match < 1.0-tol:
            bad=-1
            tgrain=tgrain+tstep
   # quit if the grain is hot enough to melt & still not hot
   # enough to satisfy the equation
        if tgrain > 10000.0:
            tgrain=-1.0
            ok=1
   # quit iterating when we are close enough
        if (match <= 1.0+tol) and (match >= 1.0-tol):
            ok=1
   # if we overshot, reduce the step size
        if oldbad != bad:
            tstep=tstep/2.0
	# CORBIN: what is this i referencing?
    t[i]=tgrain
    return t
