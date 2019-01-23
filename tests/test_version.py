# -*- coding: utf-8 -*-
"""
Unit tests for the radio_simus.version module
"""

import unittest
import sys

import radio_simus
from framework import git


try:
    import radio_simus.version
except:
    # Skip version tests for non release builds
    pass
else:
    class VersionTest(unittest.TestCase):
        """Unit tests for the version module"""

        def test_hash(self):
            githash = git("rev-parse", "HEAD")
            self.assertEqual(githash.strip(), radio_simus.version.__githash__)

        def test_version(self):
            self.assertIsNotNone(radio_simus.version.__version__)


if __name__ == "__main__":
    unittest.main()
