import unittest
import os

import numpy as np

from .. import io
from .. import footprint

desimodel_available = True
desimodel_message = "The desimodel data set was not detected."
try:
    spam = os.environ['DESIMODEL']
except KeyError:
    desimodel_available = False

class TestFootprint(unittest.TestCase):
    
    def setUp(self):
        io.reset_cache()
            
    def test_get_tile_radec(self):
        """Test grabbing tile information by tileID.
        """
        tx = io.load_tiles()
        tilefile = list(io._tiles.keys())[0]
        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [0.0, 1.0, 2.0, 3.0]
        tiles['DEC'] = [-2.0, -1.0, 1.0, 2.0]
        tiles['IN_DESI'] = [0, 1, 1, 0]
        tiles['PROGRAM'] = 'DARK'
        io._tiles[tilefile] = tiles

        #- TILEID 1 should be filtered out as not in DESI
        with self.assertRaises(ValueError):
            ra, dec = footprint.get_tile_radec(1)

        #- But TILEID 2 should be there with correct redshift
        ra, dec, = footprint.get_tile_radec(2)
        self.assertEqual((ra, dec), (1.0, -1.0))

    def test_is_point_in_desi_mock(self):
        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [0.0, 1.0, 2.0, 3.0]
        tiles['DEC'] = [-2.0, -1.0, 1.0, 2.0]
        tiles['IN_DESI'] = [0, 1, 1, 0]
        tiles['PROGRAM'] = 'DARK'

        ret = footprint.is_point_in_desi(tiles, 0.0, -2.0, return_tile_index=True)
        self.assertEqual(ret, (True, 0))

        ret = footprint.is_point_in_desi(tiles, (0.0,), (-2.0,), radius=1.605, return_tile_index=True)
        self.assertEqual(ret, ([True], [0]))

        ret = footprint.is_point_in_desi(tiles, 0.0, -3.7, radius=1.605, return_tile_index=True)
        self.assertEqual(ret, (False, 0))

        ret = footprint.is_point_in_desi(tiles, -3.0, -2.0, radius=1.605, return_tile_index=True)
        self.assertEqual(ret, (False, 0))

        ret = footprint.is_point_in_desi(tiles, tiles['RA'], tiles['DEC'], radius=1.605)
        self.assertEqual(len(ret), len(tiles))

    def test_embed_sphere(self):
        rng = np.random.RandomState(1234)
        ra = rng.uniform(0, 360, size=1000)
        dec = rng.uniform(-90, 90, size=1000)

        ret = footprint._embed_sphere(ra, dec)
        self.assertEqual(ret.shape, (1000, 3))
        norm = np.einsum('ij, ij->i', ret, ret)
        self.assertTrue(all(abs(norm - 1.0) < 1e-5))

    def test_find_points_in_tiles(self):
        rng = np.random.RandomState(1234)

        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [0.0, 1.0, 2.0, 3.0]
        tiles['DEC'] = [-2.0, -1.0, 1.0, 2.0]
        tiles['IN_DESI'] = [0, 1, 1, 0]
        tiles['PROGRAM'] = 'DARK'

        ra = rng.uniform(-10, 10, 100000)
        dec = np.degrees(np.arcsin(rng.uniform(-0.1, 0.1, 100000)))
        lists = footprint.find_points_in_tiles(tiles, ra, dec, radius=1.6058)

        # assert we've found roughly same number of objects per tile
        counts = np.array([len(i) for i in lists])
        self.assertLess(counts.std() / counts.mean(), 1e-2)

        # assert each list are indeed in the cell.
        for i, ii in enumerate(lists):
            xyz = footprint._embed_sphere(ra[ii], dec[ii])
            xyzc = footprint._embed_sphere(tiles['RA'][i], tiles['DEC'][i])
            diff = xyz - xyzc
            dist = np.einsum('ij, ij->i', diff, diff) ** 0.5
            self.assertLess(dist.max(), 2 * np.sin(np.radians(1.6058) * 0.5))

        # tiles overlapped, so we must have duplicates
        full = np.concatenate(lists)
        self.assertLess(len(np.unique(full)), len(full))

        list1 = footprint.find_points_in_tiles(tiles[0], ra, dec, radius=1.6058)
        self.assertEqual(sorted(list1), sorted(lists[0]))

    def test_find_tiles_over_point(self):
        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [0.0, 1.0, 2.0, 3.0]
        tiles['DEC'] = [-2.0, -1.0, 1.0, 2.0]
        tiles['IN_DESI'] = [0, 1, 1, 0]
        tiles['PROGRAM'] = 'DARK'

        ret = footprint.find_tiles_over_point(tiles, 0.0, -2.0)
        self.assertEqual(ret, [0, 1])

        ret = footprint.find_tiles_over_point(tiles, 1.0, -1.0, radius=1.605)
        self.assertEqual(ret, [0, 1])

        # outside
        ret = footprint.find_tiles_over_point(tiles, 0.0, -3.7, radius=1.605)
        self.assertEqual(ret, [])

        # array input
        ret = footprint.find_tiles_over_point(tiles, (0.0,), (-3.7,), radius=1.605)
        self.assertEqual(len(ret), 1)
        self.assertEqual(ret[0], [])
    
    def test_find_points_radec(self):
        """Checks if the function is successfully finding points within a certain radius of a telra and teldec"""
        import numpy as np
        answer = footprint.find_points_radec(0, 0, np.array([1.5, 0, 1.9, 0]), np.array([0, 1.5, 0, 1.3]))
        self.assertEqual(answer, [0, 1, 3])

    def test_partial_pixels(self):
        """check weights assigned to HEALPixels that partially overlap tiles"""
        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        #ADM I found a full (180) partial (406) and empty (1000) HEALPixel at nside=256
        fullpix256 = np.array([180])
        partpix256 = np.array([406])
        emptypix256 = np.array([1000])
        #ADM In the nested scheme you can traverse from 256 to 64 by integer division
        #ADM by (256/64)**2 = 16
        fullpix64 = fullpix256//16
        partpix64 = partpix256//16
        emptypix64 = emptypix256//16

        #ADM I then used:
        #ADM full = footprint.pix2tiles(256,[180])
        #ADM part = footprint.pix2tiles(256,[406])
        #ADM empty = footprint.pix2tiles(256,[20])
        #ADM to determine the tile coordinates for these pixels 
        #ADM (the first tile in  the array is full and the last tile is empty)
        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [43.05, 47.41, 47.08, 45.38]
        tiles['DEC'] = [1.54, 3.12, 3.17, 0.0]
        tiles['IN_DESI'] = [1, 1, 1, 1]
        tiles['PROGRAM'] = 'DARK'
        
        #ADM The approximate radius of DESI tiles
        radius = 1.605

        #ADM determine the weight array for pixels at nsides of 64 and 256
        pixweight64 = footprint.pixweight(64,tiles=tiles,radius=radius,precision=0.04)
        pixweight256 = footprint.pixweight(256,tiles=tiles,radius=radius,precision=0.04)
        
        #ADM check that the full pixel is assigned a weight of 1 at each nside
        self.assertTrue(np.all(pixweight64[fullpix64]==1))
        self.assertTrue(np.all(pixweight256[fullpix256]==1))

        #ADM check that the empty pixel is assigned a weight of 0 at each nside
        self.assertTrue(np.all(pixweight64[emptypix64]==0))
        self.assertTrue(np.all(pixweight256[emptypix256]==0))

        #ADM check weights of partial pixels agree reasonably at different nsides
        hirespixels = partpix64*16+np.arange(16)
        hiresweight = np.mean(pixweight256[hirespixels])
        loresweight = pixweight64[partpix64]
        #ADM really they should agree to much better than 11%. As "precision" is not set to be
        #ADM very high, this is just to check for catastrophic differences
        #ADM I checked that at precision = 0.04 this doesn't fail after 10000 attempts
        self.assertTrue(np.all(np.abs(hiresweight-loresweight) < 0.11))

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_spatial_real_tiles(self):
        tiles = io.load_tiles()
        rng = np.random.RandomState(1234)
        ra = rng.uniform(0, 360, 1000)
        dec = rng.uniform(-90, 90, 1000)

        # some points must be in the sky,
        # https://github.com/desihub/desimodel/pull/37#issuecomment-270435938
        indesi = footprint.is_point_in_desi(tiles, ra, dec, radius=1.605)
        self.assertTrue(np.any(indesi))

        # now assert the consistency between find_tiles_over_point and is_point_in_desi
        ret = footprint.find_tiles_over_point(tiles, ra, dec, radius=1.605)
        self.assertEqual(len(ret), len(ra))
        indesi2 = [len(i) > 0 for i in ret]

        # FIXME: is there a better assertion function for arrays?
        self.assertEqual(list(indesi), indesi2)

        # Just interesting to see how many tiles overlap a random point?
        ### print(np.bincount([len(i) for i in ret]))
                
def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
