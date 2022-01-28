import numpy
pro spectrumplot 

#lambda=[8.005,8.116,8.227,8.338,8.449,8.56,8.671,8.782,8.893,9.004,9.115,9.226,9.337,9.448,9.559,9.670,9.781,9.892,10.003,10.114,10.225,10.336,10.447,10.558,10.669,10.78,10.891,11.002,11.113,11.224,11.335,11.446,11.557,11.668,11.779,11.89,12.001,12.112,12.223,12.334,12.445,12.556,12.6485]
lambda_in=[2.1,3.1,3.8,4.8,8.005,8.116,8.227,8.338,8.449,8.56,8.671,8.782,8.893,9.004,9.115,9.226,9.337,9.448,9.559,9.670,9.781,9.892,10.003,10.114,10.225,10.336,10.447,10.558,10.669,10.78,10.891,11.002,11.113,11.224,11.335,11.446,11.557,11.668,11.779,11.89,12.001,12.112,12.223,12.334,12.445,12.556,12.6485,13,15,17,19,21,23,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,101]

lambdaerr=[2.1,3.1,3.8,4.8,8,8.7,8.8,9.5,9.8,10.3,10.8,11.5,11.8,11.9,12,12.5,13,19,25,60,100]

err=[0.1,0.05,0,0,0.1,0.1,0,0.05,0.1,0.1,0.05,0.05,0.2,0.05,0.1,6,0.02,0.8,5,7,8]
flux=[0.001,0.011,0.01,0.03,0.7,0.75,0.9,1.1,1.25,1.2,1.25,1.25,1.1,0.93,2.12,0.75,0.95,3.18,8.73,19.8,11.25]

x=0
fluxtot=numpy.zeroes(69)
fluxtot1=numpy.zeroes(69)
#fluxtot=numpy.zeroes(43)
#fluxtot1=numpy.zeroes(43)
#fluxtot=numpy.zeroes(21)
#fluxtot1=numpy.zeroes(21)

#old size =44

r=shift(dist(6),3,3)

#Flux without star
repeat begin
zodipic, fnu, 60, lambda_in(x), pixnum=128, nodisplay=1, dustsize=1, $
    radin=1, radout=10, starname='Beta_Pictoris', inc=80, zodis=125,alpha=0
nu=0.63*r/lambda(x)
J_one=beselj(nu,1)
PSF=(2*J_one/nu)**2
PSF[3,3]=1
k=convol(fnu,PSF)
fluxtemp=total(k)
#fluxtemp=total(fnu[120:127,*])
fluxtot(x)=fluxtemp
x=x+1
endrep until lambda_in(x) == 100

fnu=0

y=0
#Flux without star
repeat begin
zodipic, fnu, 60, lambda(y), pixnum=128, nodisplay=1, dustsize=10, $
    radin=10, radout=1000, starname='Beta_Pictoris', inc=80, zodis=3000
nu=0.63*r/lambda(y)
J_one=beselj(nu,1)
PSF=(2*J_one/nu)**2
PSF[3,3]=1
m=convol(fnu,PSF)
fluxtemp1=total(m)
#fluxtemp1=total(fnu[120:127,*])
fluxtot1(y)=fluxtemp1
y=y+1
endrep until lambda(y) == 100


fluxtot3=fluxtot1+fluxtot
#'Fluxtot_2comp1_0a.dat'  spectra0
#'Fluxtot_2comp1_0b.dat'  spectra1
#'Fluxtot_2comp1_0c.dat'  spectra2
#'Fluxtot_2comp1_0d.dat'  spectra3
#'Fluxtot_2comp1_0e.dat'  spectra4
#'Fluxtot_2comp1_0f.dat'  spectra5
#'Fluxtot_2comp1_0g.dat'  spectra6
#'Fluxtot_2comp1_0h.dat'  spectra7
#'Fluxtot_2comp1_0i.dat'  spectra8
#'Fluxtot_2comp1_0j.dat'  spectra9

save, filename='Fluxtot_2comp1_0.dat', fluxtot3, fluxtot1, fluxtot, lambda, flux

#fits_read,'lws0050_52_53_nod_stack_calib5_sub_bin3.fits',tmp,hdr
#number_of_spectra = (sxpar(hdr, 'NAXIS2')-1)/2
#wavelengths=tmp[*,0]
#spectra=tmp[*,1:number_of_spectra]
#uncertainties=tmp[*,number_of_spectra+1:number_of_spectra*2]

#normalization =66.4923, found by taking total(flux) divided by 
# avg total(spectra)

#set_plot,'ps'
#device, filename='BetaPicComp6.ps'

#window,0,color=2,title='Plots'
plot_oo, lambda, fluxtot3, xtitle='Wavelength (!7l!3m)', xstyle=1,$
    ytitle='Total Flux (Jy)',linestyle=0, xr=[1,110],yr=[1e-4,100],$
    title='Beta Pic Spectra:Inner Dust size 1!7l!3m, Outer=100!7l!3m'
#ploterror,wavelengths,spectra[*,9],uncertainties[*,9],$
#xtitle='Wavelength (!7l!3m)', ytitle='Flux (Jy)', title='Beta Pic Spectra 3.6 arcsec'

oploterr, lambdaerr, flux, err

#device,/close
#set_plot,'x'

end