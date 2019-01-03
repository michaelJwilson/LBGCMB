## Caballero 2008;  https://arxiv.org/pdf/0805.4480.pdf
import numpy                 as np
import pylab                 as pl
import pandas                as pd
import astropy.units         as u
import astropy.coordinates   as coord

from   sfd                   import ebv, reddening


c          = coord.SkyCoord(ra=[154.12, 11.1]*u.degree, dec=[-21.63, 31.65]*u.degree)

EBV        = ebv(c)
reddening  = reddening(c, survey='PS1', filters='gri')

## Work out the number of L, M and T dwarf stars in a given field.
## UDF: coordinates:    3:32:38.5,     -27:47:0.0
##      in degrees    = 53.160417,     -27.783333
##      in galactic l = 223.54060, b = -54.395430

tnum       = 2.0
field      = 'COSMOS'

surveyarea = 3.046174e-4 ## 1 degree squared in radians squared, yes!

fields     = {'COSMOS': {'b': 42.121033, 'l': 236.81763}, 'UDS': {'b':-59.964324, 'l': 169.87789}}

## Parameters of the thick disk; pc units; f_thick for reduced number density locally. 
thick_disk = {'Z_sun': 27.0, 'h_R': 2250.0, 'h_z': 500.0, 'f': 1.00}
thin_disk  = {'Z_sun': 27.0, 'h_R': 2800.0, 'h_z': 800.0, 'f': 0.04}

fname      = 'cab_table3b.lis'

# Spec-type       M_I             I-J     J-Ks   n_sol [10^-3. / pc^3]  d^b_I.
table      = pd.read_table(fname, skiprows=1, sep=r"\s*", names=['type', 'M_I', 'I-J', 'J-Ks', 'nsol', 'db'], engine='python')

print table

table['nsol'] /= 1000.0

## Scale the local density by a factor f_thinck if examining the thick disk.
table['nsol'] *= thin_disk['f']

## Convert to AB (Frei & Gunn 1995);  eh, mixing app. and abs. mag?
table['M_I'] = table['M_I'] - 0.309
table['M_J'] = table['M_I'] - table['I-J'] + 0.94
'''
## Now convert these magnitudes into z-band magnitudes
restore, 'dwarfStar_mags.sav', /verbose
Jmag  = starMags[*, 6]
zmag  = starMags[*, 13]
zminJ = -2.5*alog10(zmag/Jmag)

zminJ[23] = 3.0
M_z = M_J + zminJ

for si = 0, n_elements(zminJ) - 1 do begin

print, si, spectraltype[si], ' ', starType[si], M_J[si], M_z[si], zminJ[si]

endfor

## Calculate the minimum and maximum distances where these stars
## are a problem, so calculate the distances between which they would
## have an apparent magnitude between say 20.0 and 26.0; assume further
## selection box will be applied.
app_m_min   = 20.0
app_m_max   = 26.8

absoluteMag =  M_z

mindist = 10*10^((app_m_min - absoluteMag)/5.0) ;; pc
maxdist = 10*10^((app_m_max - absoluteMag)/5.0) ;; pc

for iii = 0, n_elements(mindist)-1 do begin
print, spectraltype[iii], " M = ", strtrim(string(absoluteMag[iii], format = '(F10.2)'), 2),  "  The min and max distances are ", string(mindist[iii], format = '(F10.2)'), string(maxdist[iii], format = '(F10.2)'), " pc"
endfor

;; quickly check that this is valid for this field:
;; zsun needs to be much smaller than dsinb

if (b ge 0) then dB = 1.0/(-(cos(b*!pi/180.0)*cos(l*!pi/180.0))/hR + (sin(b*!pi/180.0)/hz))
if (b lt 0) then dB = 1.0/(-(cos(b*!pi/180.0)*cos(l*!pi/180.0))/hR - (sin(b*!pi/180.0)/hz))

print, "The dB constant is ", dB

;; now to do some integrals.  Do an integral per spectral type:
deltamag = 0.1
n_steps = 1+ceil((app_m_max - app_m_min)/deltamag)
number_array = fltarr(n_steps, n_elements(spectraltype))
for t = 0, n_elements(spectraltype)-1 do begin
   ;; need to do discrete steps of observed apparent magnitude:
   
   app_m_array = findgen(n_steps)*deltamag + app_m_min

   for i = 0, n_steps - 1 do begin
      print, "the apparent mag range is ", app_m_array[i]- (deltamag/2.0), app_m_array[i]+ (deltamag/2.0)
      ;; find the min and max distances this corresponds to
      d1 = 10*10^((app_m_array[i] - (deltamag/2.0) - absoluteMag[t])/5.0)
      d2 = 10*10^((app_m_array[i] + (deltamag/2.0) - absoluteMag[t])/5.0)
      print, "which gives the distance range of ", d1, d2
      
      ;; check that the distance is much greater than zsun
      determinant = d1*sin(b*!pi/180.0)
      print, determinant, " pc ", d1, sin(b*!pi/180.0), b
      if (b ge 0 ) then print, "The above should be greater than Zsun = 27pc " else print, "The above should be lower than Zsun = 27pc "

      if (b ge 0) then nA = solardensity[t]*exp(-Zsun/hz)
      if (b lt 0) then nA = solardensity[t]*exp(+Zsun/hz)
      
      number = surveyarea*nA*(-((dB*d2*d2 + 2.0*dB*dB*d2+2.0*dB*dB*dB)*exp(-d2/dB)) + ((dB*d1*d1 + 2.0*dB*dB*d1+2.0*dB*dB*dB)*exp(-d1/dB)))
      print, "The number of ", spectraltype[t], " dwarfs is " , number 
      number_array[i, t] = number
   endfor
   
endfor

;; now do a plot of the overall counts per unit magnitude:

counts  = total(number_array, 2)
Mcounts = total(number_array[*, 0:6], 2)
Lcounts = total(number_array[*, 7:16], 2)
Tcounts = total(number_array[*, 17:24], 2)

;; now this is for one square degree, modify for the udf
;; udf area = 4.5 sq arcmins
;; so 1deg = 3600 sq arcmins, divide through 4.5/3600
factor = 4.5/3600.0
factor = 1.0

;; save the results to file
save, app_m_array, counts, Mcounts, Lcounts, Tcounts, number_array, spectraltype, filename = filesave
print, "Results saved to ", filesave
set_plot, 'ps'
device, filename = "plots/Dwarfstarcounts.ps", /color
plot, app_m_array, counts*factor, xthick = tnum, ythick = tnum, thick = tnum, charthick = tnum, xtitle = "z band AB mag", ytitle = "Number in 1 sq degree area per 0.1 mag", title ="Number of dwarf star contaminents in the UDS", xrange = [18,27.0], xstyle = 1
oplot,   app_m_array, Mcounts*factor, color = 70, thick = tnum
oplot,   app_m_array, Lcounts*factor, color = 150, thick = tnum
oplot,   app_m_array, Tcounts*factor, color = 200, thick = tnum

;for m = 0, 6 do begin
;oplot,  app_m_J_array, number_array[*, m]
;endfor

;;al_legend, ["M", "L", "T", "z = 7 LF", "My candidates"], linestyle = [0,0,0, 2, 0], color = [70, 150, 200, 0, 0], /bottom, /left, box = 0, thick = tnum, charthick = tnum
al_legend, ["M", "L", "T", "z = 6 LF"], linestyle = [0,0,0, 2], color = [70, 150, 200, 0], /top, /left, box = 0, thick = tnum, charthick = tnum


;restore, 'lf_8.sav'
;oplot, appm, n*factor, linestyle = 2, thick = tnum
;print, appm, n*factor

;;;;;;;;;;;;;;;;; now save results to a file for ross ;;;;;;;;;;;;;;;;;;;;
;; I need to print out magnitude, number for all dwarfs
;openw, lun, fileForRoss, /get_lun;

;printf, lun, "# J_mag  totalNumberUDF   NumberM   NumberL   NumberT"
;for i = 0, n_elements(app_m_array)-1 do begin

 ;  printf, lun, app_m_array[i], counts[i]*factor, Mcounts[i]*factor, Lcounts[i]*factor, Tcoun;ts[i]*factor

;endfor

;close, lun


;; read in my galaxy candidates:

;filecandidates = "/disk1/raab/ultravista/dualmode/catalogues/UVISTA_YJ_Y_J_sz2_iy85stack_prov9CC.cat"
;; first plot my data, read in fluxes to compute the errors easily
;readcol, filecandidates, id, ra, dec, u85, g85, r85, i85, y85, iy85stack, z85, subz2, Y, J, H, Ks, YJ, x, y, format = "L, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F"

;plothist, J, bin = 0.1, /overplot, thick = tnum

device, /close

end
'''
