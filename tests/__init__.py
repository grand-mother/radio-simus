# -*- coding: utf-8 -*-
"""
Unit tests for the radio_simus package
"""

from pathlib import Path
from os.path import split, join, realpath
import sys
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
from radio_simus import load_config

path = Path(__file__).parent / "../examples/test.config"
load_config(path)
