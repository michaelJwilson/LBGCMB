# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.focalplane.
"""
import unittest
import numpy as np
from ..focalplane import FocalPlane, generate_random_centroid_offsets, xy2radec, get_radius_mm, get_radius_deg, radec2xy
from ..focalplane import fiber_area_arcsec2, on_gfa, on_tile_gfa, get_gfa_targets
from .. import io
from astropy.table import Table

class TestFocalplane(unittest.TestCase):
    """Test desimodel.focalplane.
    """

    def test_random_offsets(self):
        """Test generating random positioner offsets.
        """
        dx, dy = generate_random_centroid_offsets(1.0)
        self.assertEqual(dx.shape, dy.shape)
        dr = np.sqrt(dx**2 + dy**2)
        self.assertAlmostEqual(np.sqrt(np.average(dr**2)), 1.0)
        self.assertLess(np.max(dr), 7)

    def test_check_radec(self):
        """Test the RA, Dec bounds checking.
        """
        F = FocalPlane()
        with self.assertRaises(ValueError):
            F._check_radec(365.0, 0.0)
        with self.assertRaises(ValueError):
            F._check_radec(0.0, 100.0)

    def test_set_tele_pointing(self):
        """Test setting the RA, Dec by hand.
        """
        F = FocalPlane()
        self.assertEqual(F.ra, 0.0)
        self.assertEqual(F.dec, 0.0)
        F.set_tele_pointing(180.0, 45.0)
        self.assertEqual(F.ra, 180.0)
        self.assertEqual(F.dec, 45.0)
    
    def test_get_radius(self):
        """Tests converting x, y coordinates on the focal plane to a radius in degrees and mm
        """
        degree = get_radius_deg(333.738, 217.766)
        truedegree = 1.5731300326614939
        radius = get_radius_mm(1.5731300326614939)
        trueradius = 398.5010456698936
        testradii = [1.5731300326614939, 1.5]
        radius1 = get_radius_mm(testradii)
        trueradius1 = [398.50104567, 378.53678987]
        self.assertAlmostEqual(truedegree, degree, 5)
        self.assertAlmostEqual(trueradius, radius, 5)
        self.assertAlmostEqual(all(trueradius1), all(radius1), 5)
    
    def test_xy2radec_new(self):
        """Tests the consistency between the conversion functions
        radec2xy and xy2radec. Also tests the accuracy of the xy2radec
        on particular cases.
        """
        truera = 8.927313423598427
        truedec = -9.324956250231294
        newra, newdec = xy2radec(8.37, -10.65, -138.345, -333.179)
        self.assertAlmostEqual(truera, newra, 5)
        self.assertAlmostEqual(truedec, newdec, 5)
        x, y = radec2xy(8.37, -10.65, newra, newdec)
        self.assertAlmostEqual(x, -138.345, 5)
        self.assertAlmostEqual(y, -333.179, 5)

    def test_xy2radec(self):
        """Test the consistency between the conversion functions
        radec2xy and xy2radec.
        """

        ra_list = np.array([0.0, 39.0, 300.0, 350.0, 359.9999, 20.0])
        dec_list = np.array([0.0, 45.0, 89.9999, -89.9999, 0.0, 89.9999])

        F = FocalPlane(ra_list, dec_list)

        # First test to check that it places x, y at the center.
        ra_obj = ra_list.copy()
        dec_obj = dec_list.copy()

        x_obj, y_obj = F.radec2xy(ra_obj, dec_obj)
        self.assertFalse(np.any(np.fabs(x_obj) > 1E-6) |
                         np.any(np.fabs(y_obj) > 1E-6),
                         ("Test Failed to recover position center with 1E-6 " +
                          "precision."))

        # Second test to check that it recovers input ra_obj,dec_obj.
        ra_out, dec_out = F.xy2radec(x_obj, y_obj)
        self.assertFalse(np.any(np.fabs(ra_out-ra_obj) > 1E-6) |
                         np.any(np.fabs(dec_out-dec_obj) >1E-6),
                         ("Test Failed to recover the input RA, Dec with " +
                          "1E-6 precision"))

    def test_area_arcsec2(self):
        """ Test area calculation
        """
        x = [0., 50., 100., 200]
        y = [-0., -50., -100., -200]
        area = fiber_area_arcsec2(x, y)
        self.assertFalse(np.any(np.fabs(area-np.array([1.97314482,  1.96207294,  1.92970506,  1.81091119])) >1E-6))

    def test_on_gfa(self):
        """Tests if on_gfa returns two lists as it is supposed to"""
        import numpy as np
        targetindices, gfalist = on_gfa(0, 0, np.array([1.0, 0, 0, 1.0]), np.array([0, 1.0, 1.5, .5]))
        self.assertEqual(0, len(targetindices))
        self.assertEqual(0, len(gfalist))
        self.assertEqual(len(targetindices), len(gfalist))

    def test_on_tile_gfa(self):
        """Tests if on_tile_gfa returns two lists as it is supposed to"""
        #- Generate some targets near a tile; some should land on the GFAs
        tile = io.load_tiles()[0]
        targets = Table()
        ntargets = 1000
        ra = np.random.uniform(tile['RA']-1.61, tile['RA']+1.61, size=ntargets)
        dec = np.random.uniform(tile['DEC']-1.61, tile['DEC']+1.61, size=ntargets)
        targets['RA'] = (ra+360) % 360
        targets['DEC'] = dec
        targetindices, gfaid = on_tile_gfa(tile['TILEID'], targets, 120)
        self.assertEqual(targetindices.size, gfaid.size)
        self.assertGreater(targetindices.size, 0)
        self.assertLess(np.max(targetindices), ntargets)
        self.assertLess(np.max(gfaid), 10)
        self.assertGreaterEqual(np.min(gfaid), 0)

        #- Test higher level get_gfa_targets
        targets['FLUX_R'] = 2000
        gfatargets = get_gfa_targets(targets, rfluxlim=1000)
        self.assertGreater(len(gfatargets), 0)
        self.assertLess(np.max(gfatargets['GFA_LOC']), 10)
        self.assertGreaterEqual(np.min(gfatargets['GFA_LOC']), 0)
        self.assertTrue(tile['TILEID'] in gfatargets['TILEID'])

        #- Shift +60 deg and verify that none are found for that tile
        targets['RA'] = (targets['RA'] + 60) % 360
        targetindices, gfaid = on_tile_gfa(tile['TILEID'], targets, 120)
        self.assertEqual(targetindices.size, 0)
        self.assertEqual(gfaid.size, 0)


def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
