========================================
Components to the throughput calculation
========================================

Introduction
============

The DESI throughput model comes from the systems engineering throughput
budget spreadsheet DESI-0347, augmented with higher resolution throughput
data for the spectrographs + CCDs from DESI-0334.  These are combined with
the KPNO extinction model ZenithExtinction-KPNO.fits and pre-calculated
fiber input geometric loss in $DESIMODEL/data/throughput/fiberloss*.dat .

The various contributions to the throughput are stored in a binary table
following the throughput format used by Specter:

https://github.com/sbailey/specter/blob/dev/doc/datamodel/throughput.md

It is basically a binary table with columns:

* wavelength : in Angstroms
* extinction : atmospheric extinction in magnitudes per airmass
* fiberinput : geometrical loss at fiber input
* throughput : all other throughput terms, e.g. mirrors, spectrograph, CCDs

Different types of sources are affected by different combinations of
throughput terms.

========== ====== === =====
Term       OBJECT SKY CALIB
========== ====== === =====
EXTINCTION yes    yes no
FIBERINPUT yes    no  no
THROUGHPUT yes    yes yes
========== ====== === =====

Source types:

* OBJECT: astronomical objects, affected by all sources of throughput loss
* SKY: sky spectra do not have a geometrical loss term for the fiber
  input.  Positioner misalignments and changes to the atmospheric
  PSF still get the same amount of sky light down the fiber.
* CALIB: calibration lamps internal to the dome do not see atmospheric
  extinction or fiber geometric loss terms.

The sections below detail the input data used to generate this table.

Atmospheric Extinction
======================

Affects astronomical objects and sky spectra, but not calibration exposures.

Depends upon airmass; extinction curve from ZenithExtinction-KPNO.fits
included in this product.

Fiber Input
===========

Affects astronomical objects, but not sky or calibration spectra.

Includes both PSF/aperture losses and pointing/guiding mis-alignment.

From:

* DESI-0347v2

  - row 16 "PSF and Aperture efficiency"
  - row 27 "Lateral errors"
  - row 52 "Fiber Defocus overall budget"
  - Multiplied then interpolated with a cubic spline

All other telescope, fiber, and instrument throughputs
======================================================

Affects all object types.

From DESI-0347v11 row 112 (total throughput) divided by the low resolution
spectrograph throughtput row 93, then multiplied by the high resolution
spectrograph+CCD throughputs in DESI-0334 *-thru*.txt files.

TO DO
=====

The spectrograph throughput numbers in DESI-0334 have been superseded by
as-built measurements from a variety of DocDB entries for each component of
the spectrographs, as listed in DESI-0347 rows 97-110.  These newer numbers
are not yet included in desimodel.
