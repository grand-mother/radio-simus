# -*- coding: utf-8 -*-
"""
Unit tests for the radio_simus.version module
"""

import unittest
import sys
import numpy as np

import radio_simus.signal_treatment as st
from framework import git


class SignalTest(unittest.TestCase):
    """Unit tests for the version module"""

    def test_add_noise(self):
        n = 100000
        vrms = 20.
        input = np.zeros((4, n))
        res = st.add_noise(vrms, input)
        for i in range(3):
            sigma = np.std(res[i, :])
            mu = np.mean(res[i, :])

        self.assertLessEqual(np.abs(mu), 5. * vrms / np.sqrt(n))
        self.assertLessEqual(np.abs(vrms) - sigma, 5.e-2 * vrms)
        self.assertEqual(res.shape, input.shape)

    def test_digitization(self):
        n = 1000

        def subtest(step, tsampling):
            input = np.zeros((4, n))
            input[0, :] = np.arange(0, n * step, step)
            input[1:, :] = np.random.standard_normal(size=(3, n))
            return st.digitization(input, tsampling), input

        step = np.random.randint(1, 10)
        tsampling = np.random.randint(1, 10) * step
        res, input = subtest(step, tsampling)
        self.assertLessEqual(res[0, 1] - res[0, 0] - tsampling, 1.e-9)
        self.assertEqual(input.shape[0], res.shape[0])
        with self.assertRaises(ValueError) :
            subtest(step, np.pi)


    def test_filter(self):
        n = 1000
        step = 1.e-9
        input = np.zeros((4, n))
        input[0, :] = np.arange(0, (n - step) * step, step)
        input[1:, :] = np.random.standard_normal(size=(3, n))
        res = st.filter(input)
        self.assertEqual(input.shape, res.shape)


if __name__ == "__main__":
    unittest.main()
