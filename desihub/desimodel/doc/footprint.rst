=========
Footprint
=========

Files
=====

* desi-tiles.par is a Yanny parameter (ASCII) file with the columns
  * tileid [1-28810] - tile number on the sky including those outside the DESI footprint
  * ra [deg] - right ascension of tile center
  * dec [deg] - declination of tile center
  * pass [1-5] - pass number, where each pass is a near-uniform tiling of the sky
  * in_desi [0-1] - 1 if in the DESI footprint, 0 otherwise
  * ebv_med [mag] - median SFD98 E(B-V) reddening value for points within the tile
  * airmass - airmass computed if observed at HA=15 deg from Kitt Peak
  * exposefac - exposure factor relative to observing the sky at zenith and E(B-V)=0
* desi-tiles.fits contains a binary table with the same information
* desi-tiles-pacman.par is the tile centers for the smaller 6-petal focal plane
* desi-tiles-pacman.fits contains a binary table with the same information
