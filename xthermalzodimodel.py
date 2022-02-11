#!/usr/bin/python
import numpy

'''
Description of code here.
'''
def xthermalzodimodel(lstar, tstar, num, inu, stepau, inc, lambda_in, radin, radout, Qabsuser, emit, lambdaQabs, iras, useralpha, userdelta, scube, userdustsize):

   # fundamental constants
    pi=3.141592
    kb=1.38066e-16
    cc=2.9979e10
    hp=6.62608e-27
    sigma=5.6705e-5
    l10=2.30259
   
   #DIRBE
   # Smooth Cloud
    alpha = 1.34
    beta = 4.14
    gamma = 0.942
    delta=0.467
    T0=286.0
    mu=0.189
    n0=1.13e-7
    em=1
   
    if (userdustsize == None):
        if (userdustsize > 3):
            lambda0=userdustsize*10
            em=(lambda0/lambda_in)**2 < 1.0
        else:
            em=Qabsuser(0)
    else:
        em=(280.0/(280.0+lambda_in))+((2.2-0.45*lambda_in) > 0)
   
    if (iras == None) :
	  # Old J Good
        alpha=1.803
        beta=4.973
        gamma=1.265
        T0=266.20
        delta=0.359
        n0=2.1527e-7  #   =1.439e-20*1.495979e13
   
    if (useralpha == None):
        alpha=useralpha
    if (userdelta == None):
        delta=userdelta
    if (scube != None):
        scube=0
   
   #************************************************************
   # More Definitions
   
   #l0 is the wavelength, in cm
    l0= lambda_in*1e-4
   
   # some math we can do now to make the Bnu calculation faster 
    chk=2.0*(kb**3.0)/((cc*hp)**2.0)
    hnuok=hp*(cc/l0)/kb
   
   #**********************************************************
   # Make some geometrical calulations 
    print('Doing the geometry')
    c1=numpy.cos(numpy.pi*inc/180.0)
    c2=numpy.sin(numpy.pi*inc/180.0)
    num2=num*2   # num2 correspondes to upixnum in zodipic
    gum=num-1
    inu=numpy.zeroes(num2, num2)
    raus=numpy.arange(num*9.0)*stepau/5.0 > 1e-8 # 9 is roughly 5 times sqrt(2)
   
    if (userdustsize is None) :
        print('Calculating the equilibrium temperature of the dust')
        #????temperaturecalc, lstar, tstar, userdustsize, raus, t, lambdaQabs, emit
    else:
        t = T0 * raus**(-delta)* (lstar**(delta/2.0))
   
   #Calculate Planck spectrum Bnu (erg s**-1 cm**-2 ster**-1 Hz**-1)
    sphereterms=hnuok/t
    sphereterms=(t**3.0)*chk*sphereterms**3.0/(numpy.exp(sphereterms)-1.0)
    t=0
	
   # multiply by a spherically symmetric term from the number density
    sphereterms=sphereterms*raus**(-alpha)
   
   # get rid of the dust outside of radout & inside of radin
    stin=(radin*5.0/stepau) > 0.0
    stout=(radout*5.0/stepau) < (num*9.0)-1
    sphereterms[0:stin]=0
    sphereterms[stout:(num*9.0)-1]=0
   
   # Fill arrays with x & z values
    xarray=numpy.zeroes(num2, num, /nozero) #???
    zarray=xarray
   #for i in range(0, num2):
    xarray[i,*]=i+0.5-num
   
    for i in range(0, num):
        zarray[*,i]=i+0.5-num
        zetacloud=abs(-c2*xarray + c1*zarray)
        xsquaredpluszsquared=xarray**2.0+zarray**2.0
   
   # Now do the azimuthally symmetric part
   #print, 'Calculating azimuthally symmetric terms.'
    zetas=numpy.arange(1001)/1000.0   # zeta goes from 0 to 1 in steps of 0.001
    if (iras is None) :
        azimuthterms=numpy.exp(-beta*(zetas**gamma))  # use for IRAS
    else:
   # DIRBE
        g= abs(zetas)-(mu/2.0)
        smallz=(zetas < mu).nonzero()	#
        smallz = smallz[0]
        g[smallz]=(zetas(smallz)**2.0)/(2.0*mu)
        azimuthterms=numpy.exp(-beta*(g**gamma))
   
    g=0
    smallz=0
    zetas=0
   
    print('Countdown')
   #*****************************************************************
    for i in range(0, gum+1):
        count=int((gum-i)/10.0)
        if (gum-i)/10.0 == count:
            print(-count)
   
	  # fill a rectangular matrix (cloud) with the distance to the center
	  # eventually cloud will contain the emission of the model at a given x 
	    y=i+0.5-num
        cloud=numpy.sqrt(xsquaredpluszsquared+y*y) > 1e-8

	  # Then make an array with zeta using a rotated z axis
        zeta=zetacloud/cloud
	  
	  # update cloud so it contains the emission, not just some geometry
        cloud=sphereterms(cloud*5.0)
        cloud=cloud*azimuthterms(zeta*1000.0)
	  
	  # clear out center of cube for iteration, if necessary
        if scube != 0 & i > gum-scube :
            gscube=scube-1
        cloud[gum-gscube:num+gscube, gum-gscube:gum]=0.0
	  
	  # Now integrate the emission along lines of sight
        inu[*,i]=stepau*sum(cloud,2)   # total along the z axis
	  # then reflect the answer to fill out inu
        inu[*,num:num2-1]=reverse(inu(*,0:gum),2)
        inu=inu+reverse(inu)
	  # Multiply the final answer by the emissivity & the local dust n<sigma>
        inu=inu*em*n0
    return
