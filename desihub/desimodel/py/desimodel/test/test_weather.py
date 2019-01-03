# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.weather.
"""
from __future__ import print_function, division

import unittest
import numpy as np
from desimodel.weather import *


class TestWeather(unittest.TestCase):
    """Test desimodel.weather.
    """
    def test_whiten_identity(self):
        """Check that whitener for samples from a Gaussian is linear.
        """
        mu, sigma = 1.2, 3.4
        gen = np.random.RandomState(seed=123)
        x = mu + sigma * gen.normal(size=100000)
        assert np.allclose(mu, np.mean(x), atol=1e-2)
        assert np.allclose(sigma, np.std(x), atol=1e-2)
        F, G = whiten_transforms(x, data_min=-1e6, data_max=+1e6)
        y = F(x)
        assert np.allclose(0, np.mean(y), atol=1e-2)
        assert np.allclose(1, np.std(y), atol=1e-2)
        ypred = (x - mu) / sigma
        assert np.allclose(ypred, y, atol=0.2)

    def test_whiten_unwhiten_roundtrip(self):
        """Check that F(G(x)) = x for a uniform distribution.
        """
        gen = np.random.RandomState(seed=123)
        lo, hi = -1.2, +3.4
        x = gen.uniform(size=10000)
        F, G = whiten_transforms(x, data_min=lo, data_max=hi)
        assert np.allclose(G(F(x)), x)

    def test_whiten_monotonic(self):
        """Check G(x) is monotonically increasing for a uniform distribution.
        """
        gen = np.random.RandomState(seed=123)
        lo, hi = -1.2, +3.4
        x = gen.uniform(size=10000)
        F, G = whiten_transforms(x, data_min=lo, data_max=hi)
        xs = np.sort(x)
        ys = F(xs)
        assert np.all(np.diff(ys) >= 0)

    def test_seeing_pdf_norm(self):
        """Check that seeing PDF is normalized.
        """
        fwhm, pdf = get_seeing_pdf()
        norm = np.sum(pdf * np.gradient(fwhm))
        assert np.allclose(norm, 1)

    def test_sample_timeseries(self):
        """Test sampling white noise with uniform 1D PDF
        """
        x = np.linspace(-1, 1, 500)
        pdf = np.ones_like(x)
        psd = lambda freq: np.ones_like(freq)
        n, nb = 1000000, 10
        gen = np.random.RandomState(1)
        xs = sample_timeseries(x, pdf, psd, n, gen=gen)
        bins, _ = np.histogram(xs, range=(-1,1), bins=nb)
        pred = n / float(nb)
        assert np.allclose(bins, pred, atol=5*np.sqrt(pred))

    def test_same_seed(self):
        """Same seed should give same samples"""
        x_grid = np.linspace(-1, 1, 500)
        pdf_grid = np.ones_like(x_grid)
        psd = lambda freq: np.ones_like(freq)
        n_sample = 1000
        gen1 = np.random.RandomState(seed=123)
        x1 = sample_timeseries(x_grid, pdf_grid, psd, n_sample, gen=gen1)
        gen2 = np.random.RandomState(seed=123)
        x2 = sample_timeseries(x_grid, pdf_grid, psd, n_sample, gen=gen2)
        self.assertTrue(np.all(x1 == x2))

    def test_different_seed(self):
        """Different seeds should give different samples"""
        x_grid = np.linspace(-1, 1, 500)
        pdf_grid = np.ones_like(x_grid)
        psd = lambda freq: np.ones_like(freq)
        n_sample = 1000
        gen1 = np.random.RandomState(seed=1)
        x1 = sample_timeseries(x_grid, pdf_grid, psd, n_sample, gen=gen1)
        gen2 = np.random.RandomState(seed=2)
        x2 = sample_timeseries(x_grid, pdf_grid, psd, n_sample, gen=gen2)
        self.assertTrue(not np.any(x1 == x2))

    def test_seeing_median(self):
        """Check that seeing has expected median
        """
        n = 100000
        gen = np.random.RandomState(seed=123)
        for m in (0.9, 1.0, 1.1, 1.2):
            x = sample_seeing(n, median_seeing=m)
            assert np.fabs(np.median(x) - m) < 0.01

    def test_transp_range(self):
        """Check that transparency has expected range
        """
        n = 100000
        gen = np.random.RandomState(seed=123)
        x = sample_transp(n)
        assert np.min(x) >= 0 and np.min(x) < 0.0001
        assert np.max(x) <= 1 and np.max(x) > 0.9999


def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
