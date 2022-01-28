#!/usr/bin/python
import numpy

'''
Description of code here.
'''


def hongphasefunction(nu):
    # output a tabulated phasefunction, pfunc, for use in
    # xscatteringzodimodel and
    # xfullzodimodel
    # pfunc is a function of cosinetheta
    costheta = (numpy.arange(1001) - 500) / 500.0
    g = [0.70, -0.20, -0.81]
    w = [0.665, 0.330, 0.005]
    nhong = 0.0795775  # 1/ ( 4 pi)
    pi = 3.1415926536
    pfunc1 = numpy.zeroes(1001)
    for i in range(0, 2 + 1):
        pfunc1 = pfunc1 + w[i] * (1.0 - g[i] * g[i]) / ((1.0 + g[i] * g[i] - 2.0 * g[i] * costheta) ** 1.5)
        pfunc1 = nhong * pfunc1
    z = numpy.zeroes(1001)
    # brightness integral for nu=1
    # as a function of cosine(theta)
    istepnum = 500.0
    for j in range(1, 999 + 1):
        eps = numpy.arccos(costheta[j])
        sineps = numpy.sin(eps)
        angi = eps + (pi - eps) * numpy.arange(istepnum + 1) / istepnum
        pfunci = numpy.zeroes(istepnum + 1)
        for i in range(0, 2 + 1):
            pfunci = pfunci + w[i] * (1.0 - g[i] * g[i]) / ((1.0 + g[i] * g[i] - 2.0 * g[i] * numpy.cos(angi)) ** 1.5)
            pfunci = nhong * pfunci
            integrand = pfunci * numpy.sin(angi)
            z[j] = (1.0 / (sineps * sineps)) * numpy.sum(integrand) * ((pi - eps) / istepnum)
    # do the conversion via equation 14 of Hong
    pfunc = pfunc1 - (nu - 1.0) * costheta * z

    return pfunc
