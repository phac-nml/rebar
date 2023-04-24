#!/usr/bin/env python3

import os

os.system(
    """
    python -m coverage run -m pytest --cov=rebar --cov-report=html --cov-report=xml \
        test/test_utils.py \
        test/test_wrappers.py
    """
)
