#!/usr/bin/python
import numpy

'''
Description of code here.
'''

def stellarparam(starname):
   
   if starname == 'Alpha_Centauri_B' or starname == 'Alf_Cen_B':
	  rstar=0.87
	  dist=1.34
	  lstar=0.45
	  tstar=5325.0 # A&A 328, 261
	  phi=0
	  gk=0.26      #Chmielewski, Y., Friel, E., Cayrel de Strobel, G., Bentolila, C.,
	  zk=0         #A&A, 263, 219
   
   if starname == 'Epsilon_Eridani' or starname == 'Eps_Eri':
	  rstar=0.98
	  lstar=0.30
	  dist=3.218 # Hipparcos
	  tstar=5156.0 # Bell, R. A., Gustafsson, B. 1989, MNRAS, 236, 653
	  phi=0
	  gk=4.75  #Drake, S. & Smith, G., ApJ, 412, 797
	  zk=-0.09 #Drake, S. & Smith, G., ApJ, 412, 797
   
   if starname == 'Alpha_Centauri_A' or starname == 'Alf_Cen_A':
	  rstar=1.23
	  dist=1.34
	  lstar=1.6
	  tstar=5770.0  # A&A 328, 261
	  phi=0
	  gk=4.3770    #unsure, 0.22 values don't work
	  zk=0
	  
   
   if starname == 'Tau_Ceti' or starname == 'Tau_Cet' :
	  tstar=5570.0
	  lstar=0.66
	  rstar=0.89
	  dist = 3.647 # Hipparcos
	  phi=0
	  gk=4.7     #Arribas, S., Martinez-Roger, C. 1988, IAU Symp.~132: The Impact
	  zk=-0.53   #of Very High S/N Spectroscopy on Stellar Physics, 132, 445
   
   if starname == 'Procyon' or starname == 'Alf_CMi' :
	  lstar=7.65
	  dist=3.497  # Hipparcos
	  tstar=6501.0 # di Benedetto, G. P. 1998, A&A, 339, 858
	  rstar=1.20 # Allen F5V
	  phi=5.51 # Mozurkewich et al. 1991
	  gk=4.00  # Fuhrmann. et al., 1997, A&A, 323, 909
	  zk=0.01  # Fuhrmann. et al., 1997, A&A, 323, 909

   
   if starname == 'Sirius' or starname == 'Alf_CMa' :
	  lstar=23.5
	  dist=2.637 # Hipparcos
	  tstar=9945.0  # di Benedetto, G. P. 1998, A&A, 339, 858
	  rstar=1.8
	  phi=5.92 # Davis & Tango 1986, Nature, 323, 234
	  gk=4.32
	  zk=0       # used sun as default
   
   if starname == 'Altair':
	  rstar=1.65
	  dist=5.143  # Hipparcos
	  tstar=8000.0
	  lstar=10.5
	  phi=0
	  gk=4.4377  # used sun as default
	  zk=0       # used sun as default
   
   if starname == 'Vega' :
	  tstar=9520.0
	  #lstar= 54.0
	  lstar= 60.0  # Backman, D.E. & Paresce, F. PP III
	  rstar=2.5
	  dist=7.756 # Hipparcos
	  phi=3.24 # Code et al. 1976, ApJ, 203, 417
	  gk=3.95  # kurget built in option
	  zk=-0.5     # unsure

   
   if starname == 'Fomalhaut' | starname == 'Alf_PsA' :
	  tstar=8800.0 # Backman, D.E. & Paresce, F. PP III
	  lstar=13.0 # Backman, D.E. & Paresce, F. PP III
	  dist=7.688 # Hipparcos
	  phi=2.10  # Code et al. 1976, ApJ, 203, 417
	  gk=3.9    # Lane, M. & Lester, J., 1984, ApJ, 281, 723
	  zk=0      # used sun default
   
   if starname == 'HR_4796' or starname == 'HR4796' :
	  rstar=1.45  # based on lstar
	  dist=67.07  # Hipparcos
	  tstar=9500.0
	  lstar=18.1 # Koerner, Ressler, Werner & Backman
	  # The big ring has an inclination of 73.1 degrees
	  # position angle 26.8 degrees see also Schneider et al. 1999
	  phi=0
	  gk=4.4377  # used sun default
	  zk=0       # used sun default
   
   if starname == 'Beta_Pictoris' or starname == 'Beta_Pic' :
	  rstar=1.46  # based on lstar
	  dist=19.3  # Hipparcos
	  tstar=8200 #Heap, S.R., Lanz, T., Hubeny, I., & Lindler, D.1994, American Astronomical Society Meeting, 26, 1389 
	  lstar=8.7 # http://astron.berkeley.edu/~kalas/disksite/pages/bpic.html
	  phi=0
	  gk=4.25   #Heap, S.R., Lanz, T., Hubeny, I., & Lindler, D.1994, American 
	  zk=-0.1   #Astronomical Society Meeting, 26, 1389 
	
   if starname == 'Fstar' :
	  rstar=1.2  # F5V Allen
	  dist=10.0  
	  tstar=6580.0 # F5V  Allen
	  lstar=2.5 
	  phi=0
	  gk=3.25
	  zk=0
   
   if starname == 'Kstar' :
	  rstar=0.7 # K7V  interpolated from Allen
	  dist=10.0 
	  tstar=3870.0 # K7V interpolated from Allen
	  lstar=0.1
	  phi=0
	  gk=4.4377
	  zk=0
   
   if starname == 'Sun' :
	  rstar=1.0
	  dist=10.0
	  tstar=5770.0
	  lstar=1.0
	  phi=0
	  gk=4.4377
	  zk=0
   
   if phi != 0:
	  print, 'Interferometrically determined angular diameter (mas):', phi
	  rsunau=6.96d10/1.49597d13
	  iradius=(phi/2000.0)*dist/rsunau
	  print, 'Combined with Hipparcos parallax gives radius (solar units):', iradius
	  rstar=iradius
   
return rstar, lstar, tstar, dist, gk, zk   
