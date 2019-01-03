====================
DESI Basis Templates
====================

Introduction
------------

This repository contains the DESI basis templates.

Versions
--------

v1.0 (2015-10-30, J. Moustakas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Original set of basis templates checked in by J. Moustakas.

Files added:

* elg_templates_v1.4.fits
* lrg_templates_v1.2.fits
* qso_templates_v1.1.fits
* star_templates_v1.2.fits
* wd_templates_v1.1.fits

v1.1 (2015-12-08, J. Moustakas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Improved wavelength binning around emission lines in the ELG templates.

Files changed:

* elg_templates_v1.5.fits

v2.0 (2016-02-22, J. Moustakas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First version of BGS templates added.  
ELG templates now utilize the high-resolution CKCz14 population synthesis models.  
Make the LRG, STAR, and WD HDU2 header FITS-compliant.  
QSO templates can now be generated on-the-fly.

Files added:

* bgs_templates_v2.0.fits
 
Files changed:

* elg_templates_v2.0.fits
* lrg_templates_v1.3.fits
* qso_templates_v2.0.fits
* star_templates_v2.0.fits
* wd_templates_v1.2.fits

v2.1 (2016-03-14, J. Moustakas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Added the first set of SNe Ia templates.

Files added:

* sne_templates_v1.0.fits

v2.2 (2016-05-10, J. Moustakas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Added a set of cooler stars and upgraded the WD templates to use the
Koester WD models.

Files changed:

* star_templates_v2.1.fits
* wd_templates_v2.0.fits

v2.3 (2017-02-27, J. Moustakas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Supplemented the DA (hydrogen atmosphere) white dwarf templates with a
set of DB (helium atmosphere) templates from Koester.  Added additional metadata
to the BGS templates so each template can be mapped to the M-XXL BGS mock
catalog. 
 
* bgs_templates_v2.1.fits
* wd_templates_v2.1.fits

v2.4 (2017-12-20, J. Moustakas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Updated LRG templates based on new target selection criteria.  See
lrg-templates.ipynb notebook in desisim/doc/nb.
 
* lrg_templates_v2.0.fits

v2.5 (2018-01-26, J. Moustakas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Updated STAR templates to v2.2.  Fixed error in wavelength array and negative
flux values (interpolation problem) in the cool stars.  
 
* star_templates_v2.2.fits
