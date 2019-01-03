# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.healpix
"""
from __future__ import print_function, division

import unittest
import numpy as np
from ..footprint import tiles2pix, tileids2pix, pix2tiles, radec2pix

try:
    import healpy
    nohealpy = False
except ImportError:
    nohealpy = True

class TestHealpix(unittest.TestCase):
    """Test healpix related routines in desimodel.footprint
    """
    @unittest.skipIf(nohealpy, 'healpy not installed')
    def test_tiles2pix_interface(self):
        """Test interface options for tiles2pix"""

        #- default tile list
        import desimodel.io
        tiles = desimodel.io.load_tiles()
        pix1 = tiles2pix(16)
        pix2 = tiles2pix(16, tiles)
        self.assertGreater(len(pix1), 0)
        self.assertListEqual(list(pix1), list(pix2))

        #- custom radius option
        pix1 = tiles2pix(16, radius=1.0)
        pix2 = tiles2pix(16, tiles, radius=1.0)
        self.assertGreater(len(pix1), 0)
        self.assertListEqual(list(pix1), list(pix2))

        #- tile IDs instead of tiles table
        pix1 = tileids2pix(16, tiles['TILEID'][0])
        pix2 = tileids2pix(16, tiles['TILEID'][0:1])
        pix3 = tiles2pix(16, tiles[0:1])
        self.assertGreater(len(pix1), 0)
        self.assertListEqual(list(pix1), list(pix2))
        self.assertListEqual(list(pix1), list(pix3))

    @unittest.skipIf(nohealpy, 'healpy not installed')
    def test_tiles2pix(self):
        """Test that the correct healpix tiles are returned for a couple values of nside.
        """
        tiles = np.zeros(3, dtype=[('RA', float), ('DEC', float)])
        tiles['RA'] = [ 335.03,  333.22,  332.35]
        tiles['DEC'] = [ 19.88,  14.84,  12.32]

        pix = tiles2pix(nside=8, tiles=tiles, radius=1.6)
        self.assertTrue( np.all(pix == np.array([196,197,208,302,303])) )

        pix = tiles2pix(nside=16, tiles=tiles, radius=1.6)
        self.assertTrue( np.all(pix == np.array([785,788,789,791,832,1209,1211,1214,1215])) )

    @unittest.skipIf(nohealpy, 'healpy not installed')
    def test_radec2pix(self):
        """test radec2pix"""
        import healpy as hp
        self.assertEqual(radec2pix(32, 0, 0), hp.ang2pix(32, np.pi/2, 0, nest=True))
        self.assertEqual(radec2pix(32, 90, 0), hp.ang2pix(32, np.pi/2, np.pi/2, nest=True))
        self.assertEqual(radec2pix(32, 45, 45), hp.ang2pix(32, np.pi/4, np.pi/4, nest=True))

    @unittest.skipIf(nohealpy, 'healpy not installed')
    def test_pix2tiles(self):
        """test roundtrip of tiles2pix -> pix2tiles"""
        import desimodel.io
        alltiles = desimodel.io.load_tiles()

        #- Test roundtrip on 5 random tiles
        nside = 16
        for tileid in np.random.choice(alltiles['TILEID'], 5, replace=False):
            for pix in tileids2pix(nside, tileid):
                tiles = pix2tiles(nside, pix)
                self.assertIn(tileid, tiles['TILEID'], '{} not in pix2tiles({},{})'.format(tileid, nside, pix))

def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
