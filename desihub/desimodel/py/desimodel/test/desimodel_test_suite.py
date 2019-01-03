# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.test.desimodel_test_suite
===================================

Used to initialize the unit test framework via ``python setup.py test``.
"""
import unittest


# This is factored out separately from runtests() so that it can be used by
# python setup.py test.
def desimodel_test_suite():
    """Returns unittest.TestSuite of desimodel tests"""
    from os.path import dirname
    desimodel_dir = dirname(dirname(__file__))
    # print(desimodel_dir)
    return unittest.defaultTestLoader.discover(desimodel_dir,
        top_level_dir=dirname(desimodel_dir))


def runtests():
    """Run all tests in desimodel.test.test_*.
    """
    # Load all TestCase classes from desimodel/test/test_*.py
    tests = desimodel_test_suite()
    # Run them
    unittest.TextTestRunner(verbosity=2).run(tests)


if __name__ == "__main__":
    runtests()
