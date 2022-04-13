#!/usr/bin/python
import numpy

'''
CORBIN: This file has some problems with the fliper functions
'''

def xscatteringzodimodel(lstar, tstar, rstar, num, inu, stepau, inc, pos, lambda_in, radin, radout, pfunc500, Qabsuser, emit, lambdaQabs, albedo, iras, scatterflag, useralpha, userdelta, scube, userdustsize):

   # fundamental constants
    pi=3.141592
    kb=1.38066e-16
    cc=2.9979e10
    hp=6.62608e-27
    sigma=5.6705e-5
    l10=2.30259
   
   #COBE
   # Smooth Cloud
    alpha = 1.34
    beta = 4.14
    gamma = 0.942
    delta=0.467
    T0=286.0
    mu=0.189
   #<n sigma> in AU**-1 measured at the Earth
    n0=1.13e-7
    al=0.18  # so sky brightness near the poles
   # this is a fit to the DIRBE emissivities
    em=(280.0/(280.0+lambda_in))+((2.2-0.45*lambda_in) > 0)
   
    if (iras == None):
        alpha=1.803
        beta=4.973
        gamma=1.265
        T0=266.20
        delta=0.359
        n0=2.1527e-7  #   =1.439e-20*1.495979e13
        em=1
   
    if (userdustsize == None):
        if (userdustsize > 3):
            lambda0=userdustsize*10
            em=(lambda0/lambda_in)**2 < 1.0
        else:
            em=Qabsuser(0)
    else:
        em=(280.0/(280.0+lambda_in))+((2.2-0.45*lambda_in) > 0)
   
    if (useralpha == None):
        alpha=useralpha
    if (userdelta == None):
        delta=userdelta
    if (scube != None):
        scube=0
    if (albedo == None):
        al=albedo
   
    num2=num*2
    gum=num-1
   
   #l0 is the wavelength, in cm
    l0= lambda_in*1e-4
   # some math we can do now to make the Bnu calculation faster 
    chk=2.0*(kb**3.0)/((cc*hp)**2.0)
    hnuok=hp*(cc/l0)/kb
   # radsolar is the radius of the sum in cm
    radsolar=6.96e10
   # rstarcm is the radius of the star in cm
    rstarcm=radsolar*rstar
   # rstarau is the radius of the star in au
    rstarau=rstarcm/1.49597e13
   # Calculate Bnu for the star
    xb=hnuok/tstar
    bnus=(tstar**3.0)*chk*xb**3.0/(numpy.exp(xb)-1.0)
    starfactor=pi * bnus * rstarau*rstarau/(stepau*stepau)
   
   #**********************************************************
   # Make some geometrical calulations 
    print('Doing the geometry')
   # work only with x < num & z < num to save time
   # then reflect the answer at the end
   
    c1=numpy.cos(pi*inc/180.0)
    s1=numpy.sin(pi*inc/180.0)
    c2=numpy.cos(pi*pos/180.0)
    s2=numpy.sin(pi*pos/180.0)
   
    inu=numpy.zeroes(num2, num2)
   
   # Do the spherically symmetric part of the physics.
    print('Calculating spherically symmetric factors')
   
   # make a vector that contains bnu as a function of radius,
   # sampled every 1/5 of a radial step
    raus=numpy.arange(num*9.0)*stepau/5.0 > 1e-8 # 9 is roughly 5 times sqrt(2)
    if (userdustsize == None) :
        print('Calculating the equilibrium temperature of the dust')
        #???temperaturecalc, lstar, tstar, userdustsize, raus, t, lambdaQabs, emit
    else:
        t = T0 * raus**(-delta)* (lstar**(delta/2.0))
   #Calculate Planck spectrum Bnu (erg s**-1 cm**-2 ster**-1 Hz**-1)
    bnu=hnuok/t
    bnu=(t**3.0)*chk*bnu**3.0/(numpy.exp(bnu)-1.0)
    t=0
   
   # compute a spherically symmetric factor from the number density
    spherefactor=raus**(-alpha)
   
   # get rid of the dust outside of radout & inside of radin
    stin=(radin*5.0/stepau)-1.0 > 0.0
    stout=(radout*5.0/stepau)+1 < num*9.0-1
    spherefactor[0:stin]=0
    spherefactor[stout:(num*9.0)-1]=0
   
   
   # Fill arrays with y & z values
    yarray=numpy.zeroes(num2, num2) # PROBLEM
    zarray=yarray
    for i in range(0, num2):
        yarray[i,0:len(yarray)]=i+0.5-num
        zarray[0:len(zarray),i]=i+0.5-num
   
   # Now do the azimuthally symmetric part
   
    print('Calculating azimuthally symmetric factors.')
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
   # save memory
    print('Countdown')
   #*****************************************************************
   
    for i in range(0, gum+1):
        count=int((gum-i)/10.0)
        if (gum-i)/10.0 == count:
            print(-count)
        x=i+0.5-num
   
   # do some serious geometry
   # start calculating the cosine of the scattering angle
        costheta=zarray*zarray   # not really costheta yet
   # rsteps will be the distance to the star in steps!
        rsteps=costheta+x*x+yarray*yarray > 1e-8   # actually this is rsteps squared
   # grab the local starlight while we have rsteps squared
        starlight=starfactor/rsteps  # like I said, it's really rsteps**2
        rsteps=numpy.sqrt(rsteps) > 1e-8    #  ah there we go...now it's rsteps
        costheta=zarray/rsteps < 0.999999 # now we have the cosine of scattering angle
   
   # Then make an array of zeta using a transformed z axis.
        zeta=abs((-s1*(s2*x + c2*yarray) + c1*zarray)/rsteps)
   
        if scatterflag == 1:
	  # compute scattered light
	  # put it all together
            cloud=spherefactor(rsteps*5.0)*azimuthterms(zeta*1000.0)*(em*bnu(rsteps*5.0) + al*starlight*pfunc500(costheta*500.0))
	  # double check that dust interior to radin is gone
            places=(rsteps < radin/stepau).nonzero()	#
            places = places[0]
            if places[0] != -1: cloud[places]=0.0
        else:
	  # no scattered light
            cloud=spherefactor(rsteps*5.0)*azimuthterms(zeta*1000.0)*em*bnu(rsteps*5.0)
	  
	  # clear out center of cube for iteration, if necessary
        if scube != 0 & i > gum-scube & i < num+scube :
            gscube=scube-1
            cloud[gum-gscube:num+gscube, gum-gscube:num+gscube]=0.0

	  # Now integrate the emission along lines of sight
        inu[i,len(inu)]=stepau*sum(cloud,2)
        # ITERATE ALONG THE ENTIRE COLUMN
	  
	  # then reflect the answer over the vertical axis to fill out inu
        if scatterflag == 1:
	  # scattered light only
            for k in gum:
                for j in len(inu):
                    inu[num:num2 - 1, len(inu)] = numpy.fliplr(inu(k, j), 1)
            for k in gum:
                for j in len(inu):
                    inu[num:num2-1,len(inu)]=numpy.fliplr(inu(k,j),1)
        else:
	  # thermal light + position angle only
            for k in gum:
                for j in len(inu):
                    inu[num:num2-1,len(inu)]=numpy.rot90(inu(k,j),2)
	  # Multiply the final answer by the local dust n<sigma>
        inu=inu*n0
        inu=numpy.rot90(inu,3)  # make positionangle=0 North
   
    return