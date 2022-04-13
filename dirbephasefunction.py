import numpy
#pro dirbephasefunction, band, pfunc

# output a tabuated phasefunction, pfunc, for use in
# xscatteringzodimodel and
# xfullzodimodel
# pfunc is a function of cosinetheta
# (the cosine of the scattering angle)
# & is tabulated for
# cos(theta)=-1.0. -0.998, -0.996 . . . 1.0

costheta=(numpy.arange(1001)-500)/500.0

pi=3.1415926536
pfunc1=numpy.zeroes(1001)

cs0=numpy.zeroes(4) 
cs1=cs0
cs2=cs0
n=cs0
## 1.25 microns
cs0[1]=-0.942
cs1[1]=0.121  # backscatter!
cs2[1]=-0.165
n[1]=3.22559  # normalization from ncalc (so integral over sphere = 1)
# 2.2 microns
cs0[2]=-0.527
cs1[2]=0.187
cs2[2]=-0.598
n[2]=0.415966
# 3.5 microns
cs0[3]=-0.431
cs1[3]=0.172
cs2[3]=-0.633
n[3]=0.324797


theta=numpy.arccos(costheta)

# CORBIN: What is band and pfunc, and why are they imported. There is a pfunc500 in zodipic.py.
pfunc=n(band)*(cs0(band) + cs1(band) * theta + numpy.exp(cs2(band)*theta))