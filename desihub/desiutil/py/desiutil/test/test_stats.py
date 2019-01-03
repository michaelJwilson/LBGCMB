# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desiutil.stats.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# The line above will help with 2to3 support.
import unittest
import numpy as np
from ..stats import perc


class TestStats(unittest.TestCase):
    """Test desiutil.stats
    """

    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_perc(self):
        """Test percentile
        """
        x = np.linspace(0, np.pi, 100)
        y = np.sin(x)
        percv = perc(y)
        np.testing.assert_allclose(percv, np.array([0.24316108649289372,
                                                   0.96590623568871437]))


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
