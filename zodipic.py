#!/usr/bin/python
import numpy
import QabsCalc
import hongphasefunction
import stellarparam
import temperaturecalc
import xfullzodimodel
import xscatteringzodimodel
import xthermalzodimodel
import xusermapzodimodel


#pro zodipic, fnu, pixsize, lambda, inclination=inclination, radin=radin, $
#    radout=radout, starname=starname, albedo=albedo, $
#    isotropic=isotropic, nodisplay=nodisplay, distance=dist, $
#    addstar=addstar, zodis=zodis, ring=ring, blob=blob, $
#    earthlong=earthlong, bands=bands, nofan=nofan, $
#    noiterate=noiterate, radring=radring, pixnum=pixnum, $
#    positionangle=positionangle, iras=iras, rstar=rstar, $
#    tstar=tstar, lstar=lstar, offsetx=offsetx, offsety=offsety, $
#    offsetz=offsetz, alpha=alpha, delta=delta, dustsize=dustsize, $
#    userdustmap=userdustmap, scaletoflux = scaletoflux

# Makes an image of what the Earth's zodiacal cloud
# would look like if it were around another star
# based on a model of the solar zodiacal cloud
# as seen by COBE

# example:  make an image of the Solar system dust viewed
# edge-on by the HST NICMOS coronagraph from 10 pc, including
# the earth ring & the dust bands, & put that image in fnu
#  zodipic, fnu, 75, 1, /ring, /bands

# example:  to make an image of what Vega would look like if it had a
# solar-type zodiacal cloud, extended out to 10 AU, as seen face-
# on by MIRLIN (a mid-IR camera) at Keck, at 12 microns, and
# include the stellar flux in the center pixel
#  zodipic, fnu, 137, 5, starname='Vega', inc=90, radout=10, /addstar

# here's a model for Epsilon Eridani as in Greaves et al. 1998, ApJ, 506,L133
#zodipic, fnu, 400, 2.2, starname='Epsilon_Eridani', radin=35, radout=70, inc=25.0, zodis=67000.0

# Also try:
#zodipic, fnu, 75, 1.1, star='HR_4796', radout=100, inc=73.1, pos=26.8, radring=70.4, /nofan, /noiter, pixnum=64, ring=10000, offsety=3

# * fnu is the output image,
#   each pixel contains a flux, in Jy
#   fnu is an array with dimensions pixnum x pixnum
#   computation time goes as pixnum**3

# * pixsize is the size of a pixel, in milliarcsec
#   You might want to oversample to get a more accurate picture of the disk.
#   If your output matrix is smaller than about 100x100, the computed
#   excess will be wrong by maybe a huge amount.

# * lambda is the wavelength, in microns

#*****************Optional Parameters**********************

# * pick a star from one of the following list using starname
#   | edit the code to add new ones
#   the Sun at 10 pc  THIS IS THE DEFAULT
#   starname='Alpha_Centauri_B' (the nearest KV star)
#   starname='Epsilon_Eridani' (the nearest KV Keck can see)
#   starname='Alpha_Centauri_A' (the nearest GV star other than the sun)
#   starname='Tau_Ceti'    (the nearest GV Keck can see, other than the sun)
#   starname='Procyon'     (the nearest FV star...Keck can see it)
#   starname='Sirius'      (the nearest AV star...Keck can see it)
#   starname='Altair'
#   starname='Vega'
#   starname='Fomalhaut'
#   starname='HR_4796'
#   starname='Beta_Pictoris'

# * radin is the inner radius of the disk, in AU
#   The default is 0.0

# * radout is the outer radius of the disk, in AU
#   Our asteroid belt stops at 3.28 AU, the 2:3 resonance with Jupiter
#   That's the default setting.
#   The DIRBE model isn't a good model of our zodiacal cloud
#   past there, but if another star has a more distant asteroid belt
#   it might be a good model for that star's disk

# *  inclination is the inclination angle in degrees
#    (0 is face-on).

# * distance is the distance to the observer, in pc
#   setting starnum overrides this parameter

# * Set addstar to add the stellar flux to the center pixel.
#   Approximates stellar spectrum as a blackbody

# * Set zodis to the number of zodis you want, & I'll multiply the
#   dust emission by that number.
#   Note that for about 1000 zodis, dust destruction by mutual collisions
#   becomes important, & this model does not include this physics

# * pixnum is the desired size of the image
#   the default is 144x144 (pixnum=144)
#   pixnum must be set to a multiple of 16
#   If you choose pixnum < 112, the model will be calaulated using
#   a internal grid that is at least 112x112 anyway

# * alpha is the radial power-law for the dust <n sigma>

# * delta is the temperature power-law for the dust <n sigma>

# The following parameters, when set, add additional detail to the
# model, & increase the running time significantly.

# * dustsize is the effective size of the dust grains, in microns
#   Setting this parameter turns on a subroutine which
#   calculates the temperature of the dust grains
#   by iteratively solving the thermal equilibrium equation.  For dust sizes
#   greater than 3 microns, it solves it assuming p=q=2 (see Backman & Paresce
#   1993, in Protostars & Planets III).  For dust sizes less than 3 microns,
#   it uses the appropriate absorption coefficient, Qabs, as calculated by
#   the program Dusty (Ivezic, NenKova, & Elitzur 1999).  Dusty uses Mie Theory
#   & optical constants for "astronomical silicate" as tabulated by Draine, B.
#   & Lee, H. (1984) ApJ, 285, 89 to calculate Qabs.  It does not affect any
#   other properties of the cloud (like optical depth).

# * set positionangle to rotate the image positionangle degrees E of N

# * Set ring=1 to add the ring & wake of dust associated with the Earth
#   according to the DIRBE model.
#   Ring is the density of the ring compared to the real ring

# * earthlong is the angle that determines the location of the
#   Earth, | rather, the Earth's wake, in degrees
#   measured in the plane of the disk from the ascending node
#   The wake trails clockwise when the disk is face-on.

# * Set bands=1 to add the bands associated with major asteroid families.

# * Set nofan to remove the fan component of the zodiacal cloud
#   Useful for visualizing the bands & ring

# * offsetx, offsety, offsetz  shift the dust from the center of the frame
#   x & y are in the plane of the disk, y points to the ascending node

# * radring is the radius of the ring#  the default (from the DIRBE model)
#   is 1.03 AU.  The width & height of the ring & blob all
#   scale with radring.

# * set noiterate if you want to force zodipic not to iterate

# * set nodisplay if you don't want zodipic to display your image

# * Set userdustmap to have zodipic calculate the emission from a map
#   of the dust distribution that you supply.  The dust density map should
#   be passed to zodipic by calling the function with the keyword
#   userdustmap equal to the 3D array containing the dust
#   distribution--units are arbitrary since the map will be scaled using
#   scaletoflux.  NOTE: ZODIPIC calculates values for userdustmap in place,
#   so it is a good idea to save userdustmap before running ZODIPIC!
#   example:  IDL>map=intarr(500,500,500)
#   IDL> zodipic, fnu, pixsize, lambda, userdustmap=map, /nodisplay
#   If radin is specified along with userdustmap, the map will be
#   cleared within this radius.  If radin is not set, the map will
#   still be cleared within the dust sublimation radius.
#   radout is ignored.  *does not yet handle scattered light*
#
#   scaletoflux is only (and must be) used in conjunction with userdustmap.
#   Set scaletoflux equal to a 2-element array containing a flux value (in
#   Jy) & a wavelength(in microns), eg
#   zodipic,....,userdustmap=map, scaletoflux=[1.34,850]
#   This tells zodipic to scale the total emission from the map to be
#   consistent with a total flux of scaletoflux[0] at a wavelength of
#   scaletoflux[1].

# Version 1 written 2/99 by
# Marc J. Kuchner
# Exoplanets & Stellar Astrophysics Laboratory
# Code 667
# NASA/Goddard Space Flight Center
# Greenbelt, MD 20771
# Marc.Kuchner@nasa.gov
# http://eud.gsfc.nasa.gov/Marc.Kuchner/home.html
# Phone: (301)286-5165  FAX: (301)286-1752

# revised 08/02 by Joannah Metz (jmetz@cfa.harvard.edu) & Sean Moran
# temperaturecalc added 3/02
# userdustmap added 6/02 by Sean Moran (smm@astro.caltech.edu)

#*************Stellar Parameters**********************
# lstar luminosity of the star in solar luminosities
# tstar is effective temperature of star, kelvin
# rstar is the radius of the star in solar radii
# dist is in pc










def zodipic():
   starname = ''
   if starname is not None:
      starname = 'Sun'

   rstar, lstar, tstar, dist, gk, zk = stellarparam(starname)
   if dist is None: udist=dist
   if rstar is None: urstar=rstar
   if tstar is None: utstar=tstar
   if lstar is None: ulstar=lstar

   # The stellar parameters we're actually going to USE are
   # ulstar, udist, urstar & utstar
   print, 'Stellar Parameters'
   print, 'Luminosity (Solar Lumin.):', ulstar, '  Distance (pc):', udist
   print, 'Radius (Solar radii):', urstar, '  Temperature (K):', utstar
   print, ' '

   #***************Fundamantal Constants**********************
   pi=3.141592
   k=1.38066d-16
   c=2.9979d10
   h=6.62608d-27
   sigma=5.6705d-5
   l10=2.30259

   # tsolar is the effective temperature of the sun 
   tsolar=5770.0
   
   # radsolar is the radius of the sum in cm
   radsolar=6.96d10
   
   # lsolar is the solar luminosity in ergs/sec
   lsolar=3.86d33
   
   #*********************More Definitions********************
   # rstarcm is the radius of the star in cm
   rstarcm=radsolar*urstar
   
   # rstarau is the radius of the star in au
   rstarau=rstarcm/1.49597d13

# wavelength in microns, of all the various DIRBE bands
# 0.5, 1.25, 2.2, 3.5, 4.9, 12.0, 25.0, 60.0, 100.0, 140.0, 240.0
   if (radout is not None): radout=3.28
   if (ring is not None): ring=0
   if (blob is not None): blob=0
   if (bands is not None): bands=0
   if (radin is not None): radin=0.0
   if (offsetx is not None): offsetx=0
   if (offsety is not None): offsety=0
   if (offsetz is not None): offsetz=0
   if (positionangle is not None): positionangle=0
   if (earthlong is not None): earthlong=0

   if (userdustmap is not None): 
       userdustmap = 0 
       scaletoflux = 0
   else:
   # if we are using a dust density map, then iterating is useless
      noiterate = 1
      zodis = 1.0 #want to make sure we dont mess up the scaling


   if (userdustmap is None) and (len(scaletoflux) != 2):
      print, ' '
      print, 'To use the userdustmap feature, please set scaletoflux, a two element vector'
      print, 'containing the desired total flux from the cloud & a fiducial wavelength.'
      print, 'I.e. scaletoflux=[desired cloud flux in Jy, wavelength in microns]'
      print, ' '

   useralpha=0
   if (alpha is None):
       useralpha=alpha
       print, 'Alpha=', useralpha
       
   userdelta=0
   userdustsize=0
   if (dustsize is None):
       userdustsize=dustsize
       print, 'Effective size of dust grains=', userdustsize
    else:
       if (delta is None):
           userdelta=delta
           print, 'Delta=', userdelta
           
# assume no oversampling
   oversample=1.0  

# upixnum is the number of pixels we're actually going to use
   if (pixnum is None):
       if pixnum/16.0 != fix(pixnum/16.0):
           print, ' '
           print, 'Please set pixnum to a multiple of 16.'
           print, ' '
   # do not allow user to use fewer than 112 pixels in the model
   # even if the output is to be smaller than that
       if pixnum < 112.0: 
          oversample=fix(112.0/pixnum)+1
       upixnum=pixnum*oversample
    else:
       # default grid size
       pixnum=144
       upixnum=144

#if using a dustmap, then adjust its size to match upixnum
   if (userdustmap is None):
       tot = total(userdustmap)
       mapsize = size(userdustmap, /dimensions)
       userdustmap = congrid(userdustmap, upixnum, upixnum, upixnum)##NNED TO FIXXXXXXXXX
       #making sure the total number density remains the same
       userdustmap = userdustmap*tot/total(userdustmap)
       pixsize = pixsize*mapsize[0]/pixnum

   print, 'The internal matrix will be', fix(upixnum), ' by', fix(upixnum)
   print, 'The output matrix will be', pixnum, ' by', pixnum
   
   upixnum=float(upixnum)
   num=upixnum/2.0
   gum=num-1
   pixsize=float(pixsize)
   radin=float(radin)
   radout=float(radout)
   stepau=pixsize*udist/(1000.0*oversample)

   if 1.415*num*stepau < radin :
       print, ' '
       print, '*************************************************'
       print,"    You are looking into the disk's central hole."
       print, '*************************************************'
       print, ' '
   
   
   if num*stepau > 3.0*radout :
       print, ' '
       print, '*************************************************'
       print,'    The disk is much smaller than the frame.'
       print, '*************************************************'
       print, ' '

   #l0 is the wavelength, in cm
   lambda=float(lambda)
   l0= lambda*1e-4
   print, 'Wavelength in microns: ', lambda

   # add scattered light if the wavelength is less than this value (in microns)
   scatterwavelength=4.2 # microns
   #*****************************************************
   # Find the radius where the dust sublimates
   delta=0.467
   T0=286.0
   if (iras is None) :
       T0=266.20
       delta=0.359
   
   tsublime=1500.0
   lambdaQabs=numpy.zeroes(112)
   a=numpy.zeroes(21)
   Qabs=numpy.zeroes(21,112)
   restore, 'QabsVars.dat'#################FIX

   if (userdustsize is None) :
       print, 'Finding the radius the dust sublimates.'
       # start at a little beyond two times the sublimation radius of a blackbody
       logsubrau=numpy.log10(2.1*(ulstar**0.5)*(T0/tsublime)**2.0)
       rok=0
       # work our way inwards till the dust temperature (tgrain)
       # exceeds the sublimation temperature
   q=2 # emission
   emit=numpy.zeroes(112)

# efficiency of absorption
# If the dust size is greater than 3 microns, the absorption is approximated
# as a power law, if it is less than 3 microns it gets the appropriate
# absorption coefficient for each wavelength from QabsCalc.

   if (userdustsize > 3.):
       # effective dust size in microns * factor
       lambda0=userdustsize*10
       emit=(lambda0/lambdaQabs)**q < 1.0
   else:
      for x in range(0, 111):
         l=lambdaQabs[x]
         Qabsreturn = QabsCalc(userdustsize, l, Qabs, lambdaQabs, a)
         emit[x]=Qabsreturn
      Qabsreturn = QabsCalc(userdustsize, lambda, Qabs, lambdaQabs, a)
      Qabsuser=Qabsreturn

   while rok == 0: # loop through radii
      tgrain = temperaturecalc(ulstar, utstar, userdustsize, 10.0**logsubrau+numpy.zeroes(1),lambdaQabs, emit)
      tgrain=tgrain[0]  # since rau & t in temperaturecalc are matrices 
      if tgrain > tsublime:
         rok=1
      else:
         logsubrau=logsubrau-0.008  # no sublimation?  decrease the trial radius by 1.859%
   rsublime=10.0**(logsubrau)
   else:
    rsublime=(ulstar**0.5)*(T0/tsublime)**(1.0/delta)
#*****************************************************
   if rsublime > radin :
      print, 'Dust sublimation temperature:', tsublime, ' K'
      print, 'Disk inner radius set to', rsublime, ' AU, .nonzero()	#
      radin=rsublime
   else:
      print, 'Inner radius, in AU:', radin
   print, 'Outer radius, in AU:', radout

   if (zodis is not None): 
      zodis=1.0
   zodis=float(zodis)
   print, 'Synthesizing an image of a Solar type zodiacal cloud x', zodis
   # Worry about whether the dust will be destroyed by mutual collisons
   alpha = 1.34
   n0=1.13d-7
   ep=1.5-alpha
   if (userdustsize is None):
       dustsizecm=dustsize*1d-4 
   else:
       dustsizecm = (5.0e-4.0) # 5 microns in diameter
   dustdensity=3.5 # grams per cc
   dustpar=3.55d-8/(dustsizecm*dustdensity)
#  Assume dust traverses the scale height of the disk twice per revolution,
# a is the heliocentric distance, & tau is the optical depth traversed
# by a dust particle since it was created at radout.
# dtau/da = dtau/dt / da/dt
# dtau/dt = 2tau0/P (P is the period=a**3/2, tau0 is the face-on optical depth) 
# da/dt = -2 dustpar/a  (see Wyatt, S. P., & Whipple, F. L. 1950, ApJ, 111, 134)  # dtau/da = -tau0 a**(-1/2) / dustpar
# rcoll is the radius where a particle released at radout
# has traversed an optical depth of unity
# The integral of dtau/da from a=radout to a=rcoll = 1
# assume tau0=0.1*n0 * rau**(1-alpha)
# & do the integral & solve for rcoll
   freepathpar=dustpar*ep/(0.1*n0*zodis)
   if freepathpar <= radout**ep :
      rcoll=(radout**ep-freepathpar)**(1.0/ep)
      print, ' '
      print, 'Heliocentric distance',1.0d4*dustsizecm
           'micron radius dust particle will be detroyed on its way in by a collision with another grain:', rcoll, ' AU'
         if radin < rcoll: 
              print, 'You might want to truncate the disk at this inner radius.'
   else:
      rcoll=0.0
      
   print, 'Step size=', stepau, ' AU'
   
   if (inclination is None):
      print, 'Inclination:', inclination, ' degrees from face-on'
      inclination=float(inclination)
   else:
      inclination=0.0
   
   if (nofan is None) :
       if ring == 0 & bands == 0 :
           print, 'You selected nofan & turned off the bands & the rings# there will be no dust.'
           zodis=0 # make sure there's no dust
        
#************************************************************

   # scatterflag makes sure you run the scattering calculation when you need it
   scatterflag=0
   if lambda < scatterwavelength :
      scatterflag=1
   if (isotropic is None) :
       pfunc500=numpy.zeroes(500)+1.0/(4.0*pi)
    else:
       print, 'Calculating phase function.'
   # calculate the phase function for our zodiacal cloud
   # use the information in Hong, S.S. 1985, A&A, 146, 67 & the fact that
   # our zodi has alpha near 1.35
      pfunc500 = hongphasefunction(1.35)
   # use some of these lines for the dirbe phase function
   # dirbeband=1 # (1 is 1.25 microns, 2 is 2.2, 3 is 3.5)
   # dirbephasefunction, dirbeband, pfunc500
   # albedo=0.204

#*********************start the action*********************
#get Fnu for the dust

# how much to magnify on each iteration:
   iterfactor=8.0

# decide how many iterations to do
# make sure there are at least minpix pixels across radin
   minpix=14.0
# minpix had better be less than upixnum!
# after one iteration there are radin/stepau steps
# after maxi iterations, there are iterfactor**maxi times as many
   maxi = fix(numpy.log10(minpix*stepau/radin)/numpy.log10(iterfactor))+1.0

# don't iterate if the user asks you not to
   if (noiterate is None): 
      maxi=1

   for i in range(1, maxi+1):
# Start with the smallest radii & move to
# larger scales with each iteration.
# how much we are magnifying by 
      magfactor=iterfactor**(maxi-i)
      print, ' '
      print, 'Iteration ', i, ' out of ', fix(maxi)
      print, 'magnification x', fix(magfactor)

# first shrink whatever we had from last iteration 
# to our new larger scale
      if i > 1 :
         fnuold=rebin(fnu,upixnum/iterfactor,upixnum/iterfactor)  # i couldn't resist the pun

# then go get a new image 
# if we are on the smallest iteration, don't leave a hole
# in the center
      scube=upixnum/(2.0*iterfactor)
      if i == 1: scube=0

#shortcut allows you to use faster models for special position angles
      shortcut=0
      if ((positionangle mod 90) == 0): 
         shortcut=1
#----------ZEROTH CASE
#this is tailored to handle the case where you supply zodipic with the
#dust distribution yourself, in the form of a 3D histogram of dust
#grain number density, stored in the variable 'userdustmap' 
      fnu = numpy.zeroes(upixnum, upixnum)
      if(userdustmap is None):
         xusermapzodimodel, ulstar, utstar, urstar, num, fnu, stepau/magfactor, \
         lambda, radin, pfunc500, userdustmap, Qabsuser, iras = iras, \
         scatterflag = scatterflag, useralpha = useralpha, userdelta = userdelta, \
         userdustsize = userdustsize, lambdaQabs = lambdaQabs, emit = emit, \
         scaletoflux = scaletoflux
      else:
         #----------FIRST CASE
         # use this program for the whole full-on salami with bands, rings etc.
         if (offsetx != 0) | (offsety != 0) | (offsetz != 0) | (ring != 0) | (bands != 0) | (scatterflag == 1 & shortcut == 0):
            xfullzodimodel, ulstar, utstar, urstar, num, fnu, stepau/magfactor, \
            inclination, positionangle, lambda, radin, radout, pfunc500, Qabsuser, emit, lambdaQabs, albedo=albedo, ring, blob,\
            earthlong, bands, nofan=nofan, offsetx, offsety, offsetz, iras=iras, scatterflag=scatterflag, useralpha=useralpha, \
            userdelta=userdelta, scube=scube, radring=radring, userdustsize=userdustsize
         else:
             #----------SECOND CASE
            if (scatterflag == 1 & shortcut == 1) :
               # use this program to do position angle = 0 at scattering wavelengths
               xscatteringzodimodel, ulstar, utstar, urstar, num, fnu, stepau/magfactor, inclination, 0.0, lambda, radin, radout, pfunc500, Qabsuser, emit,lambdaQabs, albedo=albedo, iras=iras, scatterflag=scatterflag, useralpha=useralpha, userdelta=userdelta, scube=scube, userdustsize=userdustsize
               if positionangle != 0 :
                # we took the shortcut, but we might still have some rotating to do
                # it's only by a multiple of 90 degrees, though
                  fnu=rot(fnu,positionangle) 
            #----------THIRD CASE
            else:
               if (scatterflag == 0 & shortcut == 0) :
                  # use this program to do position angle rotation for a thermal model
                  xscatteringzodimodel, ulstar, utstar, urstar, num, fnu, stepau/magfactor, inclination, positionangle, lambda, radin, radout, pfunc500, Qabsuser, emit, lambdaQabs, albedo=albedo, iras=iras, scatterflag=scatterflag, useralpha=useralpha, userdelta=userdelta,scube=scube, userdustsize=userdustsize
               else:
               #----------FOURTH CASE
               # fast no-frills, symmetrical, thermal-emission-only model
                  xthermalzodimodel, ulstar, utstar, num, fnu, stepau/magfactor, inclination, lambda, radin, radout, Qabsuser, emit, lambdaQabs, iras=iras, useralpha=useralpha, userdelta=userdelta, scube=scube, userdustsize=userdustsize
                  if positionangle != 0 :
                     # we took the shortcut, but we might still have some rotating to do
                     # it's only by a multiple of 90 degrees, though
                     fnu=rot(fnu,positionangle) 

# Compute total flux from this iteration, in Jy
      if (scaletoflux is not None): #since scaletoflux already in Jy
         fnuit = total(fnu)*1d23*((pixsize/oversample)/(magfactor*1000.0*206265.0))**2.0
      else: 
         fnuit = total(fnu)
      print, 'Flux added this iteration:', fnuit 
      # add the old image to the new image
      if i > 1:
         gscube=upixnum/(2.0*iterfactor)-1.0
         print, 'gscube=',gscube
      fnu[gum-gscube:num+gscube,gum-gscube:num+gscube]=fnu[gum-gscube:num+gscube,gum-gscube:num+gscube]+fnuold
   # end i iteration loop

   if (scaletoflux is not None): 
    # Convert disk surface brightness from cgs to Jy/ster
       fnu=fnu*1d23
    # now convert that to Jy
       fnu=fnu*((pixsize/oversample)/(1000.0*206265.0))**2.0
       fnu=fnu*zodis

   fdisk=total(fnu)
   print, 'Total flux (Jy):', fdisk

# if we oversampled, bin back down
   if upixnum != pixnum :
       fnu=rebin(fnu, pixnum, pixnum)
       fnu=fnu*fdisk/total(fnu)   # make sure that the total flux is correct


# redefine num because we are now working again with a 
# pixnum x pixnum grid
   num=pixnum/2.0

#Calculate flux from star from Planck spectrum 
# Bnu (erg s**-1 cm**-2 ster**-1 Hz**-1)
   nu = c/l0
   xb=h*nu/(k*utstar)
   bnu=xb**3.0/(numpy.exp(xb)-1.0)
   bnu=bnu*2.0*((k*utstar)**3.0)/((c*h)**2.0)

   distcm=udist* (3.085678e18.0) # distance to star in cm
   fstar = 1e23 * pi * bnu * ((rstarcm/distcm)**2.0)
   
   index=-1
   # Add star to the image, if required
   if (addstar is None) :
       if (2*rstarau) > stepau :
           print, "Warning: Star is bigger than one pixel."
           fnu[num-1,num-1]=fnu[num-1,num-1]+fstar/4.0
           fnu[num,num-1]=fnu[num,num-1]+fstar/4.0
           fnu[num-1,num]=fnu[num-1,num]+fstar/4.0
           fnu[num,num]=fnu[num,num]+fstar/4.0
           print, 'The Stellar flux has been divided evenly among pixels'
           print, fix(num-1), ',', fix(num-1)
           print, fix(num), ',', fix(num-1)
           print, fix(num-1), ',', fix(num)
           print, fix(num), ',', fix(num)

   print, 'Stellar flux', fstar, ' Jy'
   if (fstar != 0):
      print, 'Excess (Disk flux/Stellar Flux)', fdisk/fstar

   # Display the result FIX
   # if the grid is small, make it look bigger
   case 1 of
   num <= 9: rfactor=16 
   (num <= 18) 
   (num > 9) : rfactor=8
   #(num <= 37) 
   (num > 18) : rfactor=8
   #(num <= 73) 
   (num > 37) : rfactor=4
   #(num <= 145) 
   (num > 73) : rfactor=2
   #else: rfactor=1
   endcase

   if pixnum < 600 & (nodisplay is not None) :
      pic=rebin(fnu, pixnum*rfactor, pixnum*rfactor, /sample)
      sz=size(pic)
      pic=rotate(pic,2)  # to make North up in the display
   
   #DISPLAY
   wtitle='Surface Brightness'
   #window, 0, title=wtitle,xsize=sz(1),ysize=sz(2)
   
   wtitle='Log Surface Brightness'
   #window, 2, title=wtitle,xsize=sz(1),ysize=sz(2)
   places=(pic > 0).nonzero()	#
   places = places[0]
   if places[0] != -1 :
       amin=min(pic(places))
       amin=amin > 1e-20
       pic=pic > amin
       
return