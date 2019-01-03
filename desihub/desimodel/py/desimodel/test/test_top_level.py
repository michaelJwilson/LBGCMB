# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
test top-level desimodel functions
"""
import unittest
import re
from .. import __version__ as theVersion


class TestTopLevel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.versionre = re.compile(r'''([0-9]+!)?   # Epoch
                                       ([0-9]+)     # Major
                                       (\.[0-9]+)*  # Minor...
                                       ((a|b|rc|\.post|\.dev)[0-9]+)?''',
                                   re.X)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_version(self):
        """Ensure the version conforms to PEP386/PEP440.
        """
        try:
            self.assertRegex(theVersion, self.versionre)
        except AttributeError:
            self.assertRegexpMatches(theVersion, self.versionre)


def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
