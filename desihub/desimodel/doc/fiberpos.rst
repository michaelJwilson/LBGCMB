=================
Fiber Positioners
=================

Input Files
===========

DESI-0530 defines the location of positioners on the focal plane.

DESI-2721 defines the mapping from cassettes of 50 positioners/fibers on the
    focal plane to bundles of fibers on the spectrograph slit heads

Output Files
============

Within a cassette, the fiber order is randomized.
desimodel.inputs.fiberpos.update() randomizes the fibers within a cassette
and outputs $DESIMODEL/data/focalplane/fiberpos.[fits, txt, png] with the
mapping of positioner -> fiber number.  fiberpos-all.[fits, ecsv] also
includes non-spectrograph fiber positioner locations such as fiducials
and sky monitors.

Coordinate Systems
==================

As defined in DESI-0481v1 table 1 (last row "CS5 DESI Focal Plane"),
when the telescope is parked at zenith the +x direction is East (+RA),
+z is pointed toward the ground, and thus +y points south (-dec).

Note that y and z are sign swapped wrt to the other DESI coordinate systems:
most coordinate systems think about +z pointing toward the sky, while the
focal plane thinks about +z as pointing away from the focal plane and thus
away from the sky.
