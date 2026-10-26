"""Shared pytest fixtures / path setup for the wdm_av_builder tests."""

from __future__ import annotations

import os
import sys

# make the package importable when pytest is run from anywhere
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
