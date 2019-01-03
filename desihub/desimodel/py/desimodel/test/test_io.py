# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.io.
"""
import os
import sys
import uuid
import tempfile
import numpy as np
from astropy.table import Table
import unittest
from .. import io
#
# Try to import specter.
#
specter_available = True
specter_message = "The specter package was not detected."
try:
    import specter
except ImportError:
    specter_available = False
#
# Try to load the DESIMODEL environment variable
#
desimodel_available = True
desimodel_message = "The desimodel data set was not detected."
try:
    spam = os.environ['DESIMODEL']
except KeyError:
    desimodel_available = False
    specter_available = False
    specter_message = desimodel_message


class TestIO(unittest.TestCase):
    """Test desimodel.io.
    """
    @classmethod
    def setUpClass(cls):
        global specter_available, desimodel_available
        cls.specter_available = specter_available
        cls.desimodel_available = desimodel_available
        cls.tempdir = tempfile.mkdtemp(prefix='testio-')
        cls.trimdir = os.path.join(cls.tempdir, 'trim')
        cls.testfile = os.path.join(cls.tempdir, 'test-abc123.fits')

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.tempdir):
            import shutil
            shutil.rmtree(cls.tempdir)

    def setUp(self):
        """Ensure that any desimodel.io caches are clear before running
        any test.
        """
        io.reset_cache()
        if os.path.exists(self.testfile):
            os.remove(self.testfile)

    def tearDown(self):
        pass

    def test_reset_cache(self):
        """Test cache reset (two examples at least)
        """
        self.assertTrue(io._fiberpos is None)
        self.assertTrue(isinstance(io._tiles, dict))
        self.assertEqual(len(io._tiles), 0)
        fiberpos = io.load_fiberpos()
        tiles = io.load_tiles()
        self.assertTrue(io._fiberpos is not None)
        self.assertTrue(isinstance(io._tiles, dict))
        self.assertEqual(len(io._tiles), 1)
        io.reset_cache()
        self.assertTrue(io._fiberpos is None)
        self.assertTrue(isinstance(io._tiles, dict))
        self.assertEqual(len(io._tiles), 0)

    @unittest.skipUnless(specter_available, specter_message)
    def test_load_throughput(self):
        """Test loading of throughput files.
        """
        for channel in ('b', 'r', 'z'):
            t = io.load_throughput(channel)

    @unittest.skipUnless(specter_available, specter_message)
    def test_load_psf(self):
        """Test loading of PSF files.
        """
        for channel in ('b', 'r', 'z'):
            t = io.load_psf(channel)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_desiparams(self):
        """Test loading of basic DESI parameters.
        """
        p = io.load_desiparams()
        self.assertTrue(isinstance(p, dict))
    
    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_gfa(self):
        """Test loading of gfa location data.
            """
        gfa = io.load_gfa()
        # length of the GFA table should be 40 since there are each corner of the 10 GFAs are included
        self.assertEqual(len(gfa), 40)
        for key in ('PETAL', 'X', 'Y', 'Z', 'Q', 'RADIUS_DEG', 'RADIUS_MM'):
            self.assertIn(key, gfa.dtype.names)
        
    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_fiberpos(self):
        """Test loading of fiber positioner data.
        """
        fiberpos = io.load_fiberpos()
        self.assertEqual(len(fiberpos), 5000)
        # Check workaround for astropy 1.0.x bug for lower -> upper col names.
        for key in ('FIBER', 'LOCATION', 'SPECTRO', 'X', 'Y', 'Z'):
            self.assertIn(key, fiberpos.dtype.names)
            x = fiberpos[key]

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_platescale(self):
        """Test loading platescale.txt file.
        """
        p1 = io.load_platescale()
        p2 = io.load_platescale()
        self.assertTrue(p1 is p2)  #- caching worked

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_targets(self):
        """Test loading of tile files.
        """
        data = io.load_target_info()
        #- Test a few keys, but not everything
        self.assertIn('ntarget_lrg', data.keys())
        self.assertIn('nobs_elg', data.keys())
        self.assertIn('success_qso', data.keys())

    @unittest.skip('Skip *until* the pixel weights file is in the DESIMODEL directory')
#    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_pix_file(self):
        """Test loading of the file of HEALPixel weights.
        """
        nside = 16
        pix = io.load_pixweight(nside)
        npix = len(pix)
        #ADM Test the length of the returned array is appropriate to the passed nside
        self.assertEqual(npix,12*nside*nside)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_tiles(self):
        """Test loading of tile files.
        """
        # starting clean
        self.assertEqual(io._tiles, {})
        t0 = io.load_tiles(cache=False)
        self.assertEqual(len(io._tiles), 0)
        # loading tiles fills the cache with one items
        t1 = io.load_tiles(onlydesi=False)
        self.assertEqual(len(io._tiles), 1)
        tile_cache_id1 = id(list(io._tiles.items())[0])
        # reloading, even with a filter, shouldn't change cache
        t2 = io.load_tiles(onlydesi=True)
        self.assertEqual(len(io._tiles), 1)
        tile_cache_id2 = id(list(io._tiles.items())[0])
        self.assertEqual(tile_cache_id1, tile_cache_id2)
        #- Temporarily support OBSCONDITIONS as u2 (old) or i4 (new)
        self.assertTrue(np.issubdtype(t1['OBSCONDITIONS'].dtype, 'i4') or \
                        np.issubdtype(t1['OBSCONDITIONS'].dtype, 'u2') )
        self.assertTrue(np.issubdtype(t2['OBSCONDITIONS'].dtype, 'i4') or \
                        np.issubdtype(t2['OBSCONDITIONS'].dtype, 'u2') )
        self.assertLess(len(t2), len(t1))
        # All tiles in DESI are also in full set.
        self.assertTrue(np.all(np.in1d(t2['TILEID'], t1['TILEID'])))
        # I think this is the exact same test as above, except using set theory.
        self.assertEqual(len(set(t2.TILEID) - set(t1.TILEID)), 0)
        t3 = io.load_tiles(onlydesi=False)
        tile_cache_id3 = id(list(io._tiles.items())[0])
        self.assertEqual(tile_cache_id1, tile_cache_id3)
        self.assertTrue(np.issubdtype(t3['OBSCONDITIONS'].dtype, 'i4') or \
                        np.issubdtype(t3['OBSCONDITIONS'].dtype, 'u2') )
        # Check for extra tiles.
        a = io.load_tiles(extra=False)
        self.assertEqual(np.sum(np.char.startswith(a['PROGRAM'], 'EXTRA')), 0)
        b = io.load_tiles(extra=True)
        self.assertGreater(np.sum(np.char.startswith(b['PROGRAM'], 'EXTRA')), 0)
        self.assertLess(len(a), len(b))

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_tiles_alt(self):
        # starting clean
        self.assertEqual(io._tiles, {})
        t1 = Table(io.load_tiles())
        t1.write(self.testfile)
        #- no path; should fail since that file isn't in $DESIMODEL/data/footprint/
        if sys.version_info.major == 2:
            with self.assertRaises(IOError):
                t2 = io.load_tiles(tilesfile=os.path.basename(self.testfile))
        else:
            with self.assertRaises(FileNotFoundError):
                t2 = io.load_tiles(tilesfile=os.path.basename(self.testfile))

        #- with path, should work:
        t2 = Table(io.load_tiles(tilesfile=self.testfile))
        self.assertTrue(np.all(t1 == t2))

        #- cache should have two items now
        self.assertEqual(len(io._tiles), 2)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_tiles_consistency(self):
        """Test consistency of tile files.

        - Validate different tile loading schemes.
        - Test numerical consistency between FITS & ECSV versions.
        - Assure that PROGRAM column has no trailing whitespace.
        """
        fitstiles = io.findfile('footprint/desi-tiles.fits')
        ecsvtiles = io.findfile('footprint/desi-tiles.ecsv')
        # tf=TilesFits  tt=TilesTable  te=TilesEcsv
        tf = io.load_tiles(onlydesi=False, extra=True)
        tt = Table.read(fitstiles)
        te = Table.read(ecsvtiles, format='ascii.ecsv')
        self.assertEqual(sorted(tf.dtype.names), sorted(tt.colnames))
        self.assertEqual(sorted(tf.dtype.names), sorted(te.colnames))

        for col in tt.colnames:
            self.assertTrue(np.all(tf[col]==tt[col]), 'fits[{col}] != table[{col}]'.format(col=col))
            if np.issubdtype(tf[col].dtype, float):
                self.assertTrue(np.allclose(tf[col], te[col], atol=1e-4, rtol=1e-4), 'fits[{col}] != ecsv[{col}]'.format(col=col))
            else:
                self.assertTrue(np.all(tf[col]==te[col]), 'fits[{col}] != ecsv[{col}]'.format(col=col))

        self.assertTrue(not np.any(np.char.endswith(tf['PROGRAM'], ' ')))
        self.assertTrue(not np.any(np.char.endswith(tt['PROGRAM'], ' ')))
        self.assertTrue(not np.any(np.char.endswith(te['PROGRAM'], ' ')))

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_trim_data(self):
        '''Test trimming data files for lightweight tests'''
        psffile = io.findfile('specpsf/psf-b.fits')
        if os.path.getsize(psffile) > 20e6:
            from .. import trim
            indir = os.path.join(os.getenv('DESIMODEL'), 'data')
            trim.trim_data(indir, self.trimdir)
            self.assertTrue(os.path.isdir(self.trimdir))
            self.assertGreater(len(list(os.walk(self.trimdir))), 1)


def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
