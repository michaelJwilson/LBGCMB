# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import unittest
import os

from ..inputs import gfa

desimodel_available = True
desimodel_message = "The desimodel data set was not detected."

class TestInputs(unittest.TestCase):
    """Test desimodel.inputs
    """
    def setUp(self):
        pass
    
    def test_build_gfa_table(self):
        gfatable = gfa.build_gfa_table('testname.ecsv')
        self.assertTrue(os.path.exists('testname.ecsv'), "Test Failed to create a valid file")
        os.remove('testname.ecsv')
        
        
def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
    
