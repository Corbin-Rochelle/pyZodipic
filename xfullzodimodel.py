#!/usr/bin/python
import numpy

'''
Description of code here.
'''


def xfullzodimodel(lstar, tstar, rstar, num, inu, stepau, inc, pos, lambda_in, radin, radout, pfunc500, Qabsuser, emit, lambdaQabs, albedo, ring, blob, earthlong, bands, nofan, offsetx, offsety, offsetz, iras,scatterflag, useralpha, userdelta, scube, radring, userdustsize):

    #fundamental constants
    pi=3.1415926536
    pi2=pi*2.0
    kb=1.38066e-16
    cc=2.9979e10
    hp=6.62608e-27
    sigma=5.6705e-5
    l10=2.30259
   
    # DIRBE
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
    em=1
   
    if userdustsize == None:
        if userdustsize > 3:
            lambda0=userdustsize*10
            em=(lambda0/lambda_in)**2 < 1.0
        else:
            em=Qabsuser(0)
    else:
        em=(280.0/(280.0+lambda_in))+((2.2-0.45*lambda_in) > 0)
   
    if iras == None:
   # Old J Good
        alpha=1.803
        beta=4.973
        gamma=1.265
        T0=266.20
        delta=0.359
        n0=2.1527e-7
   
    if useralpha == None:
        alpha=useralpha
    if userdelta == None:
        delta=userdelta
    if albedo == None:
        al=albedo
    if scube != None:
        scube=0
   
   # ***************************************************************
   #  Solar Ring
    nsr=1.83e-8
    rsr=1.03
    if radring == None:
        rsr=radring
    sigrsr=0.025*rsr/1.03  # make width of ring scale with the radius
    sigrsr2=2.0*sigrsr*sigrsr
    sigzsr=0.054*rsr/1.03  # make the height of the ring scale with the radius
   # Trailing Blob
    ntb=1.9e-8
    rtb=1.06*rsr/1.03  # make the location of the blob scale with the ring
    sigrtb=0.10*rsr/1.03  # make width of blob scale with the ring
    sigrtb2=2.0*sigrtb*sigrtb
    sigztb=0.091*rsr/1.03  # make height of the blob scale with the ring
    heltb=-10.0*numpy.pi/180.0  # converted to radians
    sigheltb=12.1*numpy.pi/180.0  # converted to radians
    sigheltb2=sigheltb*sigheltb
   
   #***************************************************************
   # Dust bands
    nb1=5.59e-10
    delb1=8.78*numpy.pi/180.0
    vb1=0.10
    pb1=4.0
    rcut1=1.5
    nb2=1.99e-9
    delb2=1.99*numpy.pi/180.0
    vb2=0.90
    pb2=4.0
    rcut2=0.94
    nb3=1.44e-10
    delb3=15.0*numpy.pi/180.0
    vb3=0.05
    pb3=4.0
    rcut3=1.5
    rcut120=rcut1**20.0
    rcut220=rcut2**20.0
    rcut320=rcut3**20.0
    delb16=delb1**6.0
    delb26=delb2**6.0
    delb36=delb3**6.0
    delb14=delb1**4.0
    delb24=delb2**4.0
    delb34=delb3**4.0
   
   #************************************************************
   # More Definitions
   
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
   
   # rstarau is the radius of the star in AU
    rstarau=rstarcm/1.49597e13
   
    deltax=offsetx/stepau
    deltay=offsety/stepau
    deltaz=offsetz/stepau
   
   # make sure earthlong is between 0 & 2*pi
    earthlong=(earthlong*numpy.pi/180.0) % pi*2
    if earthlong < 0:
        earthlong = earthlong + pi2
   
   # Calculate Bnu for the star
    xb=hnuok/tstar
    bnus=(tstar**3.0)*chk*xb**3.0/(numpy.exp(xb)-1.0)
   
    starfactor=numpy.pi * bnus * rstarau*rstarau/(stepau*stepau)
   
   #**********************************************************
   # Make some geometrical calulations 
    print('Doing some geometry')
   
    c0=numpy.cos(pos*numpy.pi/180.0 )
    s0=numpy.sin(pos*numpy.pi/180.0 )
    c1=numpy.cos(inc*numpy.pi/180.0)
    s1=numpy.sin(inc*numpy.pi/180.0)
    c2=numpy.cos(0.5*numpy.pi + earthlong)
    s2=numpy.sin(0.5*numpy.pi + earthlong)
   
    inu=numpy.zeroes(num2, num2)
   
   # Do the spherically symmetric part of the physics.
    print('Calculating spherically symmetric factors')
   
   # make a vector that contains bnu as a function of radius,
   # sampled every 1/5 of a radial step
    raus=numpy.arange(num*9.0)*stepau/5.0 > 1e-8 # 9 is roughly 5 times sqrt(2)
   
    if userdustsize == None:
        print('Calculating the equilibrium temperature of the dust')
	  #temperaturecalc, lstar, tstar, userdustsize, raus, t, lambdaQabs, emit
    else:
        t = T0 * raus**(-delta)* (lstar**(delta/2.0))
   #Calculate Planck spectrum Bnu (erg s**-1 cm**-2 ster**-1 Hz**-1)
    bnu=hnuok/t
    bnu=(t**3.0)*chk*bnu**3.0/(numpy.exp(bnu)-1.0)
    t=0
	
   # compute a spherically symmetric factor from the number density
    spherefactor=raus**(-alpha)
    if nofan == None:
        spherefactor=spherefactor*0.0
   
   # Fill arrays with y & z values
    yarray=numpy.zeroes(num2, num2)
    zarray=yarray
    for i in range(0, num2):
        yarray[i,0:len(yarray)]=i+0.5-num
        zarray[0:len(zarray),i]=i+0.5-num
   
   # do some math with these matrices for the coordinate transformation
    trans1=-(s0*c1*c2 + c0*s2)*yarray + s1*c2*zarray + c2*deltax -s2* deltay
    trans2=-(s0*c1*s2 - c0*c2)*yarray + s1*s2*zarray + s2*deltax +c2*deltay
    trans3=s0*s1*yarray + c1*zarray+deltaz
   
    print('Calculating azimuthally symmetric factors.')
    zetas=numpy.arange(1001)/1000.0   # zeta goes from 0 to 1 in steps of 0.001
    if iras == None:
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
    for i in range(0, num2):
        count=int(((num2-1)-i)/10.0)
        if ((num2-1)-i)/10.0 == count:
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
   
        x3=(c0*c1*c2 - s0*s2)*x + trans1
        y3=(c0*c1*s2 + s0*c2)*x + trans2
        z3=-c0*s1*x + trans3
   
        z3=abs(z3)
        rstepst=numpy.sqrt(x3*x3+y3*y3+z3*z3)
        zeta=abs(z3/rstepst)
   
        raut=stepau*rstepst  # in units of AU
        if ring != 0 | bands != 0:
            zaut=stepau*z3 # in units of AU
   
   # make an array to put the numberdensity of the ring, wake + bands in
        nd=zarray*0.0
   
        if ring != 0 :
   # Add the Earth Ring
            nd=nd+ring*nsr*numpy.exp(-(((raut-rsr)**2.0)/(sigrsr2) + zaut/sigzsr ))
   
        if blob != 0 :
    # Add the Earth Blob
   # no inclination | delta
   # We're going to need this angle
            deltahel=numpy.atan(y3,x3)
            dhel=heltb-deltahel
            places=(dhel > numpy.pi).nonzero()	#
            places = places[0]
        if places[0] != -1: dhel[places]=dhel[places] - pi2
        nd=nd+blob*ntb*numpy.exp(-(((raut-rtb)**2.0)/(sigrtb2)+zaut/sigztb+(dhel*dhel/(2.0*sigheltb2)) ))
   
   # Add the dust bands
        if bands != 0 :
   #print, 'Adding the Dust Bands'
            zau6=zaut**6.0
            zau4=zaut**4.0
            rau20=raut**20.0
            nd=nd+3.0*bands*nb1*numpy.exp(-(zau6/delb16))*(vb1+(zau4/delb14))*(1.0-numpy.exp(-(rau20/rcut120)))
            nd=nd+3.0*bands*nb2*numpy.exp(-(zau6/delb26))*(vb2+(zau4/delb24))*(1.0-numpy.exp(-(rau20/rcut220)))
            nd=nd+3.0*bands*nb3*numpy.exp(-(zau6/delb36))*(vb3+(zau4/delb34))*(1.0-numpy.exp(-(rau20/rcut320)))
		 
        if scatterflag == 1:
   # compute scattered light
   # put it all together
            cloud=(spherefactor[rstepst*5.0]*azimuthterms(zeta*1000.0)*n0+nd)*(em*bnu(rsteps*5.0) + al*starlight*pfunc500(costheta*500.0))
        else:
   # no scattered light
            cloud=(spherefactor[rstepst*5.0]*azimuthterms(zeta*1000.0)*n0+nd)*(em*bnu(rsteps*5.0))
   
   # get rid of dust interior to radin & exterior to radout
        places=(raut < radin).nonzero()	#
        places = places[0]
        if places[0] != -1:
            cloud[places]=0.0
        places=(raut > radout).nonzero()	#
        places = places[0]
        if places[0] != -1:
            cloud[places]=0.0
   
   # clear out center of cube for iteration, if necessary
        if scube != 0 & i > gum-scube & i < num+scube :
            gscube=scube-1
            cloud[gum-gscube:num+gscube, gum-gscube:num+gscube]=0.0
		 
   # Now integrate the emission along lines of sight
        inu[i,len(inu)]=stepau*sum(cloud,2)
    return