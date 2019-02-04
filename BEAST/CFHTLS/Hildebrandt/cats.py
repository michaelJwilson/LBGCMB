import os
import pylab as pl 
import astropy.io.fits as fits


'''
Selection criteria:

--  Colour / colour.
    E.g.  'u_g_ISO_cor2' and 'g_r_ISO_cor2'


--  dat[1].data['MASK_u']      == 0
--  dat[1].data['MASK_g']      == 0
--  dat[1].data['MASK_r']      == 0
--  dat[1].data['MASK_i']      == 0
--  dat[1].data['MASK_z']      == 0
--  dat[1].data['MASK_STARS']  == 0 
--  dat[1].data['masksa']      == 0
--  dat[1].data['maskgscw2']   == 0
--  dat[1].data['CLASS_STAR'] < 0.9
'''

root     = os.environ['CSCRATCH']
dat      = fits.open(root + '/Hildebrandt/CFHTLS/D1.fits')

good   = (dat[1].data['MASK_u'] == 0) & (dat[1].data['MASK_g'] == 0) & (dat[1].data['MASK_r'] == 0) & (dat[1].data['MASK_i'] == 0) & (dat[1].data['MASK_z'] == 0)\
       & (dat[1].data['MASK_STARS'] == 0) & (dat[1].data['masksa'] == 0) & (dat[1].data['maskgscw2'] == 0) & (dat[1].data['CLASS_STAR'] < 0.9)

umg    = dat[1].data['u_g_ISO_cor2'] 
gmr    = dat[1].data['g_r_ISO_cor2']
rmi    = dat[1].data['r_i_ISO_cor2']
imz    = dat[1].data['i_z_ISO_cor2']

udrop  = (1.5 < umg) & (-1.0 < gmr) & (gmr < 1.2) & (1.5 * gmr < umg -0.75)
ucat   = good & udrop

gdrop  = (1.0 < gmr) & (-1.0 < rmi) & (rmi < 1.0) & (1.5 * rmi < gmr -0.80)
gcat   = good & gdrop

rdrop  = (1.2 < rmi) & (-1.0 < imz) & (imz < 0.7) & (1.5 * imz < rmi -1.00)
rcat   = good & rdrop

umag   = dat[1].data['MAG_ISOCOR_u'] 
gmag   = dat[1].data['MAG_ISOCOR_g']
rmag   = dat[1].data['MAG_ISOCOR_r']
imag   = dat[1].data['MAG_ISOCOR_i']
zmag   = dat[1].data['MAG_ISOCOR_z']

'''
##  u-drops.
pl.plot(gmr[ucat], umg[ucat], 'kx') 

pl.xlabel(r'$(g-r)$')
pl.ylabel(r'$(u-g)$')

pl.xlim(-0.5, 4.0)
pl.ylim(-0.5, 4.0)

pl.show()

pl.clf()


##  g-drops
pl.plot(rmi[gcat], gmr[gcat], 'kx')

pl.xlabel(r'$(r-i)$')
pl.ylabel(r'$(g-r)$')

pl.xlim(-0.5, 4.0)
pl.ylim(-0.5, 4.0)

pl.show()


##  r-drops                                                                                                                                                                                                                                                                            
pl.plot(imz[rcat], rmi[rcat], 'kx')

pl.xlabel(r'$(i-z)$')
pl.ylabel(r'$(r-i)$')

pl.xlim(-0.5, 4.0)
pl.ylim(-0.5, 4.0)

pl.show()
'''
pl.clf()

##  Photometric redshifts.                                                                                                                                                                                       
HypZCWW = dat[1].data['Z_PHOT']          ##  HyperZ run with the CWW.
HypZBC3 = dat[1].data['Z_PHOT_BC']       ##  HyperZ run with the BC03.  
BPZ     = dat[1].data['Z_B_V3']          ##  BpZ.

for cat, band in zip([ucat, gcat, rcat], ['u', 'g', 'r']): 
  pl.clf()
  pl.hist(dat[1].data['Z_PHOT'][cat],    bins=25, label='%s-drop HPZ-CWW' % band, alpha=0.4)
  pl.hist(dat[1].data['Z_PHOT_BC'][cat], bins=25, alpha=0.4, label='HPZ-BC03')
  pl.hist(dat[1].data['Z_B_V3'][cat],    bins=25, alpha=0.4, label='BPZ')

  pl.xlabel(r'$z$')
  pl.ylabel(r'$dN/dz$')

  pl.legend()
  pl.show()

print('\n\nDone.\n\n')
