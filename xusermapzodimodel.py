#!/usr/bin/python
import numpy

'''
Description of code here.
'''

def xusermapzodimodel(lstar, tstar, rstar, num, inu, stepau, lambda_in, radin, pfunc500, im, Qabsuser, iras, scatterflag, useralpha, userdelta, userdustsize, lambdaQabs, emit, scaletoflux):
#*****************************************************************
# fundamental constants
    pi=3.1415926536
    pi2=pi*2.0
    kb=1.38066e-16
    cc=2.9979e10
    hp=6.62608e-27
    sigma=5.6705e-5
    l10=2.30259
   
    alpha = 1.34
    beta = 4.14
    gamma = 0.942
    delta=0.467
    T0=286.0
    mu=0.189
   #<n sigma> in AU**-1 measured at the Earth
    n0=1.13e-7
    al=0.36  # so Hong phase function agrees with DIRBE at 90 degrees
   
    if [iras is None] :
        alpha=1.803
        beta=4.973
        gamma=1.265
        T0=266.20
        delta=0.359
        n0=2.1527e-7  #   =1.439e-20*1.495979e13

    if useralpha != 0:
        alpha=useralpha
    if userdelta != 0:
        delta=userdelta

    num2=num*2
#l0 is the wavelength, in cm
    l0= lambda_in*1e-4

# some math we can do now to make the Bnu calculation faster 
    chk=2.0*[kb**3.0]/[[cc*hp]**2.0]
    hnuok=hp*[cc/l0]/kb
    hnuok_scale = hp*[cc/[scaletoflux[1]*1.0e-4]]/kb

# If the dust size is less than 3 microns, then it takes the appropriate
# absorption coefficient & stores it to `em'
    em=1
    if (userdustsize is None):
        if (userdustsize > 3):
            lambda0=userdustsize*10
            em=[lambda0/lambda_in]**2 < 1.0
        else:
            em=Qabsuser(0)
    else:
        em=[280.0/[280.0+lambda_in]]+[[2.2-0.45*lambda_in] > 0]

# where we're gonna put the results
    inu=numpy.zeroes(num2, num2)

# make a vector that contains bnu as a function of radius,
# sampled every 1/5 of a radial step

    raus=[numpy.arange(num*9.0)+0.5]*stepau/5.0 > 1e-8
#gotta find temperatures of dust before we get Bnu
    if [userdustsize is None] :
        print('Calculating the equilibrium temperature of the dust')
	  #???temperaturecalc, lstar, tstar, userdustsize, raus, t, lambdaQabs, emit
    else:
        t = T0 * raus**(-delta)*(lstar**(delta/2.0))
#create a cubic array containing distance from center, then
#temperature, then bnu at that position
    print('generating emission map from the dust density map...')
    x = [numpy.arange(num)+0.5]*stepau
    y = make_array(num, value = 1.0)
    plane = x#y
    plane = plane**2.0+[transpose(plane)]**2.0
    plane=[[rotate(plane,2),rotate(plane,3)],[rotate(plane,1),plane]]
    cube = make_array(num2, num2, num2, value = 1.0)
    for i in range(0, num-1):
        cube[num-i-1, :, :] =plane[:, :]+[stepau*[float(i)+0.5]]**2.0
        cube[num+i, :, :] = cube[num-i-1, :, :]
    cube = sqrt(cube)

   #clear density map within radin
    indices = [cube[[num-radin-3]:[num+radin+2][num-radin-3]:[num+radin+2] [num-radin-3]:[num+radin+2]] < radin cnt].nonzero()	#
    indices = indices[0]
    #??= len(indices)
    if cnt != 0:
        im[indices] = 0.0
   
   #fill with temperatures
    cube[:, :, :]= t[floor(cube[:, :, :]/stepau*5.0)]
   
   #scale the density map to match the known flux scaletoflux[0] at the
   #wavelength scaletoflux[1]
    if len(scaletoflux) == 2 and scaletoflux[1] != lambda_in:
	  #calculate the dust emission coeff at the scaleto wavelength
        print('Scaling to match ', scaletoflux[0], 'Jy at ', scaletoflux[1], 'microns')
        em2=1.0
        if userdustsize == None:
            if userdustsize > 3:
                em2=[280.0/[280.0+scaletoflux[1]]]+[[2.2-0.45*scaletoflux[1]] > 0]
            else:
                em2=Qabsuser(0)
        else:
            em2=[280.0/[280.0+scaletoflux[1]]]+[[2.2-0.45*scaletoflux[1]] > 0]
   #calculate the total emission
        flux = hnuok_scale/cube
        flux = [cube**3.0]*chk*flux**3.0/[numpy.exp(flux)-1.0]
        scalefactor = total(flux*im*em2)
        im = im*scaletoflux[0]/scalefactor
   
   #Calculate Planck spectrum Bnu (erg s**-1 cm**-2 ster**-1 Hz**-1)
   # for grains at each of the temperatures in t_cube
    bnu=hnuok/cube
    bnu=[cube**3.0]*chk*bnu**3.0/[numpy.exp(bnu)-1.0]
   
   #multiply the density map by bnu & emission coefficient
    im = bnu*im*em
   
   #if scaling to  flux at the same wavelength as the image we're making,
   #take this shortcut.
    if scaletoflux[1] == lambda_in:
        im = im/total(im)*scaletoflux[0]
   
   # Now integrate the emission along lines of sight
    if len(scaletoflux) == 2:
        inu = total(im, 3)
    else:
        inu=stepau*total(im,3)
    return