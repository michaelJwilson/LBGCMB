# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.install.
"""
from os.path import abspath, dirname
import unittest
from .. import __version__ as desimodel_version
from ..install import default_install_dir, svn_export, install


class TestInstall(unittest.TestCase):
    """Test desimodel.install.
    """

    def test_default_install_dir(self):
        """Test setting default install directory.
        """
        d1 = default_install_dir()
        d2 = dirname(dirname(dirname(dirname(dirname(dirname(__file__))))))
        self.assertEqual(abspath(d1), abspath(d2))

    def test_svn_export(self):
        """Test svn export command.
        """
        base_url = "https://desi.lbl.gov/svn/code/desimodel/{0}/data"
        cmd = svn_export()
        self.assertEqual(cmd[2], base_url.format('trunk'))
        cmd = svn_export('trunk')
        self.assertEqual(cmd[2], base_url.format('trunk'))
        cmd = svn_export('branches/v3')
        self.assertEqual(cmd[2], base_url.format('branches/v3'))
        cmd = svn_export('1.2.3')
        self.assertEqual(cmd[2], base_url.format('tags/1.2.3'))

    def test_install(self):
        """Test the install function.
        """
        pass


def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
