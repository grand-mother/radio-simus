# -*- coding: utf-8 -*-
"""
Unit tests for the radio_simus package
"""

from pathlib import Path
from radio_simus import load_config

path = Path(__file__).parent / "../examples/test.config"
load_config(path)
