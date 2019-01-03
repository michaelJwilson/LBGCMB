"""
tests desispec.sky
"""

import unittest
import pdb

import numpy as np
import os
from desispec.frame import Frame
#from desispec.qa import QA_Frame, QA_Exposure, QA_Brick, QA_Prod
from desispec.qa.qa_frame import QA_Frame
from desispec.qa.qa_exposure import QA_Exposure
from desispec.qa.qa_brick import QA_Brick
from desispec.qa.qa_prod import QA_Prod
from desispec.io import write_qa_frame, write_qa_brick, load_qa_frame, write_qa_exposure, findfile, write_frame
from desispec.io import write_fiberflat, specprod_root
from desispec.test.util import get_frame_data, get_calib_from_frame, get_fiberflat_from_frame
#from uuid import uuid4
from shutil import rmtree

class TestQA(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.nspec = 6
        cls.nwave = 20
        id = 1
        cls.night = '20160101'
        cls.expid = 1
        # Paths
        os.environ['DESI_SPECTRO_REDUX'] = os.environ['HOME']
        os.environ['SPECPROD'] = 'desi_test_qa'
        # Files
        cls.testDir = specprod_root()
        cls.qafile_b0 = findfile('qa_data', night=cls.night, expid=cls.expid, specprod_dir=cls.testDir, camera='b0')
        cls.qafile_b1 = findfile('qa_data', night=cls.night, expid=cls.expid, specprod_dir=cls.testDir, camera='b1')
        cls.qafile_exp = cls.testDir+'/exposures/'+cls.night+'/{:08d}/qa-{:08d}'.format(id,id)
        cls.qafile_brick = cls.testDir+'/brick/3582m005/qa-3582m005.yaml'
        cls.flux_pdf = cls.testDir+'/exposures/'+cls.night+'/{:08d}/qa-flux-{:08d}.pdf'.format(id,id)
        # Files for exposure fibermap QA figure
        cls.frame_b0 = findfile('frame', night=cls.night, expid=cls.expid, specprod_dir=cls.testDir, camera='b0')
        cls.frame_b1 = findfile('frame', night=cls.night, expid=cls.expid, specprod_dir=cls.testDir, camera='b1')
        cls.fflat_b0 = findfile('fiberflat', night=cls.night, expid=cls.expid, specprod_dir=cls.testDir, camera='b0')
        cls.fflat_b1 = findfile('fiberflat', night=cls.night, expid=cls.expid, specprod_dir=cls.testDir, camera='b1')
        cls.exp_fmap_plot = cls.testDir+'/test_exp_fibermap_plot.png'

    @classmethod
    def tearDownClass(cls):
        """Cleanup in case tests crashed and left files behind"""
        for filename in [cls.qafile_b0, cls.qafile_b1, cls.flux_pdf, cls.frame_b0, cls.frame_b1,
                         cls.fflat_b0, cls.fflat_b1, cls.exp_fmap_plot]:
            if os.path.exists(filename):
                os.remove(filename)
                #testpath = os.path.normpath(os.path.dirname(filename))
                #if testpath != '.':
                #    os.removedirs(testpath)
        if os.path.exists(cls.testDir):
            rmtree(cls.testDir)

    def _make_frame(self, camera='b0', flavor='science', night=None, expid=None, nspec=3, objtype=None):
        if night is None:
            night = self.night
        if expid is None:
            expid = self.expid
        # Generate
        frame = get_frame_data(nspec=nspec, objtype=objtype)
        frame.meta = dict(CAMERA=camera, FLAVOR=flavor, NIGHT=night, EXPID=expid)
        return frame

    def _write_flat_files(self):
        # Frames
        fb0 = self._make_frame(camera='b0', flavor='flat', nspec=10, objtype='FLAT')
        _ = write_frame(self.frame_b0, fb0)
        fb1 = self._make_frame(camera='b1', flavor='flat', nspec=10, objtype='FLAT')
        _ = write_frame(self.frame_b1, fb1)
        # Fiberflats
        ff0 = get_fiberflat_from_frame(fb0)
        write_fiberflat(self.fflat_b0, ff0)
        ff1 = get_fiberflat_from_frame(fb1)
        write_fiberflat(self.fflat_b1, ff1)

    def _write_qaframes(self):
        """Write QA data frame files"""
        frm0 = self._make_frame(camera='b0')
        frm1 = self._make_frame(camera='b1')
        qafrm0 = QA_Frame(frm0)
        qafrm1 = QA_Frame(frm1)
        # SKY
        qafrm0.init_skysub()
        qafrm1.init_skysub()
        qafrm0.qa_data['SKYSUB']['METRICS'] = {}
        qafrm1.qa_data['SKYSUB']['METRICS'] = {}
        qafrm0.qa_data['SKYSUB']['METRICS']['NSKY_FIB'] = 10
        qafrm1.qa_data['SKYSUB']['METRICS']['NSKY_FIB'] = 30
        # FLUX
        qafrm0.init_fluxcalib()
        qafrm1.init_fluxcalib()
        qafrm0.qa_data['FLUXCALIB']['METRICS'] = {}
        qafrm0.qa_data['FLUXCALIB']['METRICS']['ZP'] = 24.
        qafrm0.qa_data['FLUXCALIB']['METRICS']['RMS_ZP'] = 0.05
        qafrm1.qa_data['FLUXCALIB']['METRICS'] = {}
        qafrm1.qa_data['FLUXCALIB']['METRICS']['ZP'] = 24.5
        qafrm1.qa_data['FLUXCALIB']['METRICS']['RMS_ZP'] = 0.05
        # WRITE
        write_qa_frame(self.qafile_b0, qafrm0)
        write_qa_frame(self.qafile_b1, qafrm1)

    def _write_qabrick(self):
        """Write a QA data brick file"""
        qabrck = QA_Brick()
        # ZBEST
        qabrck.init_zbest()
        qabrck.data['ZBEST']['METRICS'] = {}
        qabrck.data['ZBEST']['METRICS']['NFAIL'] = 10
        write_qa_brick(self.qafile_brick, qabrck)

    def test_init_qa_frame(self):
        #- Simple Init call
        qafrm1 = QA_Frame(self._make_frame(flavor='science'))
        assert qafrm1.flavor == 'science'

    def test_init_qa_fiberflat(self):
        #- Init FiberFlat dict
        qafrm = QA_Frame(self._make_frame(flavor='flat'))
        qafrm.init_fiberflat()
        assert qafrm.qa_data['FIBERFLAT']['PARAMS']['MAX_RMS'] > 0.

        #- ReInit FiberFlat dict
        qafrm.init_fiberflat(re_init=True)
        assert qafrm.qa_data['FIBERFLAT']['PARAMS']['MAX_RMS'] > 0.

    def test_init_qa_fluxcalib(self):
        #- Init FluxCalib dict
        qafrm = QA_Frame(self._make_frame(flavor='science'))
        qafrm.init_fluxcalib()
        assert qafrm.qa_data['FLUXCALIB']['PARAMS']['MAX_ZP_OFF'] > 0.

        #- ReInit FluxCalib dict
        qafrm.init_fluxcalib(re_init=True)
        assert qafrm.qa_data['FLUXCALIB']['PARAMS']['MAX_ZP_OFF'] > 0.

    def test_init_qa_skysub(self):
        #- Init SkySub dict
        qafrm = QA_Frame(self._make_frame(flavor='science'))
        qafrm.init_skysub()
        assert qafrm.qa_data['SKYSUB']['PARAMS']['PCHI_RESID'] > 0.

        #- ReInit SkySub dict
        qafrm.init_skysub(re_init=True)
        assert qafrm.qa_data['SKYSUB']['PARAMS']['PCHI_RESID'] > 0.

    def test_qa_frame_write_load_data(self):
        # Write
        frm0 = self._make_frame()
        qafrm0 = QA_Frame(frm0)
        write_qa_frame(self.qafile_b0, qafrm0)
        # Load
        qafrm2 = load_qa_frame(self.qafile_b0, frm0)
        assert qafrm2.night == qafrm0.night


    def test_init_qa_exposure(self):
        """Test simple init.
        """
        from os import environ
        cache_env = {'SPECPROD': None, 'DESI_SPECTRO_REDUX': None}
        for k in cache_env:
            if k in environ:
                cache_env[k] = environ[k]
            environ[k] = './'
        qaexp = QA_Exposure(1, '20150211', 'arc')
        self.assertEqual(qaexp.expid, 1)
        for k in cache_env:
            if cache_env[k] is None:
                del environ[k]
            else:
                environ[k] = cache_env[k]

    def test_qa_exposure_load_write_data(self):
        #- Test loading data
        self._write_qaframes()
        qaexp = QA_Exposure(self.expid, self.night, 'science', specprod_dir=self.testDir)
        assert 'b0' in qaexp.data['frames']
        assert 'b1' in qaexp.data['frames']
        # Write
        write_qa_exposure(self.qafile_exp, qaexp)

    def test_exposure_fibermap_plot(self):
        from desispec.qa.qa_plots import exposure_fibermap
        self._write_flat_files()
        exposure_fibermap('b', self.expid, 'meanflux', outfile=self.exp_fmap_plot)

    """
    # This needs to run as a script for the figure generation to pass Travis..
    def test_qa_exposure_fluxcalib(self):
        #- Perform fluxcalib QA on Exposure (including figure)
        self._write_qaframes()
        qaexp = QA_Exposure(1, self.night, specprod_dir=self.testDir,
                            flavor='dark')
        qaexp.fluxcalib(self.flux_pdf)
    """

    def test_init_qa_brick(self):
        #- Simple Init calls
        qabrck = QA_Brick(name='tst_brick')
        assert qabrck.brick_name == 'tst_brick'
        #
        qabrck.init_zbest()
        assert qabrck.data['ZBEST']['PARAMS']['MAX_NFAIL'] > 0

    def test_init_qa_prod(self):
        qaprod = QA_Prod(self.testDir)

    def test_qa_frame_plot(self):
        from desispec.qa import qa_plots
        from desispec.qa import qa_frame
        # Frame
        frame = get_frame_data(500)
        # Load calib
        fluxcalib = get_calib_from_frame(frame)
        # QA Frame
        tdict = {}
        tdict['20190829'] = {}
        dint = 20
        tdict['20190829'][dint] = {}
        tdict['20190829'][dint]['flavor'] = 'science'
        tdict['20190829'][dint]['b'] = {}
        tdict['20190829'][dint]['b']['FLUXCALIB'] = {}
        tdict['20190829'][dint]['b']['FLUXCALIB']['METRICS'] = {}
        tdict['20190829'][dint]['b']['FLUXCALIB']['METRICS']['BLAH'] = 1
        qaframe = qa_frame.QA_Frame(tdict)
        # Plot
        qa_plots.frame_fluxcalib('tmp.pdf', qaframe, frame, fluxcalib)

    def runTest(self):
        pass


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)



#- This runs all test* functions in any TestCase class in this file
if __name__ == '__main__':
    unittest.main()
    #qa_frame_plot_test()
