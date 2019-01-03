=======================
desimodel Release Notes
=======================

0.9.1 (unreleased)
------------------

* No changes yet

0.9.0 (2017-09-19)
------------------

* Added desimodel.focalplane.radec2xy, which converts RA, Dec coordinates to x, y coordinates on the focal plane, which accepts vector inputs.
* Added desimodel.focalplane.on_gfa() and its respective helper functions to check if a target is on a GFA of arbitrary telescope pointing
* Added desimodel.focalplane.on_tile_gfa() to check return a list of indices of targets on a specific tile
* Added desimodel.focalplane.get_gfa_targets() to return a table with added columns GFA_LOC and TILEID that consists of all targets on any GFA on any tile satisfying a minimum flux in the r-band. 
* Unittests for the desimodel.focalplane functions were updated accordingly.
* Added desimodel.footprint.find_points_in_tel_range() to return a list of indices withnin a radius of an arbitray telescope pointing, unaware of tiles (Added respective unittest)
* Adds desimodel.focalplane.fiber_area_arcsec2()
* Updates tests to work with trimmed data subset

0.8.0 (2017-08-07)
------------------

* Add new weather module to specify assumed atmospheric seeing and transparency
  distributions at KPNO, with accompanying DESI-doc and jupyter notebook.
* Remove seeing module, which is superseded by new weather module.
* Added `desimodel.footprint.pixweight()` in :mod:`desimodel.footprint` to create an array of what fraction
  of every HEALPixel at a given nside overlaps the DESI footprint
* Also added `desimodel.footprint.tiles2fracpix()` to estimate which HEALPixels overlap the footprint edges
* Added `desimodel.io.load_pixweight()` in :mod:`desimodel.io` to load the array created by
  `desimodel.footprint.pixweight()` and resample it to any HEALPix nside
* Modified path to Lya SNR spectra files used in desi_quicklya.py, used in Lya Fisher forecast.
* Added desimodel.inputs.build_gfa_table and its helper functions to write a .ecsv file for GFA data
* Added desimodel.io.load_gfa to return the GFA data table
* Added desimodel.focalplane.xy2radec, which converts x,y coordinates on the focal plane to RA, Dec coordinates
* don't print warnings in desimodel.io if specter isn't installed

0.7.0 (2017-06-15)
------------------

* Added desimodel.footprint.tiles2pix and .pix2tiles for mapping healpix
  to DESI tiles.
* fixed psf-quicksim.fits units to be astropy-friendly
* added `desimodel.io.load_target_info()`

0.6.0 (2017-03-27)
------------------

* Add desimodel.seeing module with functions that model the expected DESI
  zenith seeing at 6355A, with an accompanying jupyter notebook.
* Altered xy offset RMS calculation in focalplane.py to scale the distribution
  RMS rather than the sample standard deviation.
* Update focal plane to positioner mapping
* z-channel 250 um CCD instead of 500 um CCD
* Update DocDB -> desimodel update method for fiberpos and throughput

0.5.1 (2016-12-01)
------------------

* By default, desimodel.io.load_tiles now excludes PROGRAM=EXTRA layers
* Adds desi-tiles.* tests

0.5.0 (2016-11-21)
------------------

* Moved test of focalplane code into the actual test suite.
* Preparing for Python 3.
* Changed default svn version to trunk and added error handling to
  :command:`install_desimodel_data`.
* Update template module file to reflect DESI+Anaconda infrastructure.
* Add code to generate random centroid offsets in :mod:`desimodel.focalplane`.
* Add jupyter notebook documenting new throughput files of `PR#29`_.
* Use Astropy-recommended method of reading FITS data tables.
* Remove reference to Travis scripts in MANIFEST.in.

.. _`PR#29`: https://github.com/desihub/desimodel/pull/29

0.4.5 (2016-07-15)
------------------

* Fixed a minor bug that made the help message for :command:`install_desimodel_data`
  garbled.
* Add additional files to lightweight test data to work with quickgen

0.4.4 (2016-03-15)
------------------

* Allow :command:`desiInstall` to download and install the data from svn.
* No changes to data in svn.

0.4.3 (2016-03-10)
------------------

* "First" post-separation tag.
* Added :func:`desimodel.trim.trim_data` for trimming a data directory into a
  lightweight version for testing.
* svn data includes targets.dat: preliminary numbers for MWS and BGS densities
  (Still waiting upon supporting technote).

0.4.2 (2016-02-04)
------------------

* Improved svn download instructions in the README file.
* Changes to data on svn side

  * updated desi.yaml with dark vs. bright exptime
  * updated targets.dat to include MWS placeholders

* :func:`desimodel.io.load_desiparams` adds 'exptime' -> 'exptime_dark' key for temporary
  backwards compatibility
* Removed deprecated fibers module
* Use `ci-helpers`_ to handle most of the dirty work of Travis build scripts.
* Make `specter`_ import errors more verbose.

.. _`ci-helpers`: https://github.com/astropy/ci-helpers
.. _`specter`: https://github.com/desihub/specter

0.4.1 (2016-01-25)
------------------

* Last tag prior to separating desimodel into code (GitHub) and data (svn)
  repositories.
* pip install support (BAW).
* Replace fitsio dependency with astropy.io.fits.

0.4 (2015-12-14)
----------------

* Added tile file for the bright time survey.

0.3.8 (2015-10-30)
------------------

* Adds python io library (SJB).

0.3.7 (2015-04-16)
------------------

* Tag to support dogwood production (SJB).

0.3.6 (2015-01-30)
------------------

* Adds :meth:`desimodel.focalplane.FocalPlane.xy2radec` from Jaime (SJB).

0.3.5 (2014-12-28)
------------------

data/targets/targets.dat
    added fractions for sky and stdstar fibers (SJB).

py/desimodel/focalplane.py
    bug fixes for transformations (SJB).

0.3.4 (2014-09-23)
------------------

* Fix a simple import error (BAW).

0.3.3 (2014-09-23)
------------------

* Fix a simple version error (BAW).

0.3.2 (2014-09-23)
------------------

* Change how version is set (BAW).
* Updated target numbers.

0.3.1 (2014-07-23)
------------------

* Also updated quicksim sn-spec* file output, using IDL version which is slightly
  more optimistic than the python version (diff is dark current?) (SJB).

0.3 (2014-07-23)
----------------

* Updated throughput files for real.
* Added initial "compare_versions.py" script to make it easier to visualize
  differences in versions.  This script should grow as various parameters
  change; right now it only makes a thoughput difference plot (SJB).
* Updated throughput files from 0334v3 (spectro) and 0347v5 (system throughput)
  Correction: thoughput files didn't make it into that change (SJB, 2014-07-08).
* Updated psf-b.fits and psf-quicksim.fits to match new npix_y for blue
  STA/ITL CCDs (SJB, 2014-07-08).

0.2 (2014-07-08)
----------------

2014-07-07 SJB
~~~~~~~~~~~~~~

* Added ELG spectrum with continuum and multiple emission lines

2014-07-07 David Kirkby
~~~~~~~~~~~~~~~~~~~~~~~

Python quicksim

* add readnoise contributions in quadrature during the downsampling
* Refactor for speed, results now named ndarray, updated plots
* Allow different base directories

2014-07-02 DJS
~~~~~~~~~~~~~~

* Put sky back to dimmer UVES sky model

0.1 (2014-07-01)
----------------

2014-06-29 SJB
~~~~~~~~~~~~~~

* Extended fiberloss range from 3500-10000 instead of 3600-10000
* Added data/throughput/fiberloss-qso.dat (same as fiberloss-star.dat)

2014-06-27 SJB
~~~~~~~~~~~~~~

* Updated data/focalplane/platescale.txt with latest from DESI-0329v14.
  This includes a new "theta" column.
* Updated desi.yaml from DESI-0347v4.  This removes the FWHM and wavemin/max
  params which are not derived quantities associated with the PSFs.
* Updated throughput files with new numbers from DESI-0347v4.
* Updated spectrograph throughput files with new numbers from DESI-0334v2.
* Updated py/fiberloss.py -> bin/fiberloss.py .  Biggest change is ELG
  half light radius 0.35" -> 0.45" which drops us below 7-sigma.
* Updated data/throughput/fiberloss-\*.dat files with calculation based
  upon fiberloss.py
* bin/psf2quicksim.py extracted PSF parameters needed for quicksim.
    - pro/desi_quicksim.pro updated, but it still treats FWHM as constant
      rather than wavelength dependent.
    - python quicksim will be broken until it is updated to use new inputs.
* Reorganized data/inputs/throughput/
* spots2psf.py: leftover spot mirroring bug removed, PSFs updated

2014-06-12 SJB
~~~~~~~~~~~~~~

* Updated throughputs to not double count central obscuration.
* Updated PSF files to remove throughputs to avoid possible inconsistency.
* Added wavemin_all, wavemax_all to desi.yaml with min/max wavelength
  seen by all spectra

2014-06-06 SJB
~~~~~~~~~~~~~~

* Updated CCD pixel dimensions and regenerated PSFs to match.
