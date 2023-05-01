#!/usr/bin/env python3

import subprocess
import sys
import shutil
import os

tmp_dir = "test/tmp"

if os.path.exists(tmp_dir):
    shutil.rmtree(tmp_dir, ignore_errors=True)

cmd_str = (
    "python -m coverage run -m pytest --cov=rebar --cov-report=html --cov-report=xml"
    " test/test_wrappers.py"
    # `test_edge_cases`` depends on `test_wrappers`` to be run first
    " test/test_edge_cases.py"
    # `test_utils` is the biggest test suite, and takes the longest(?)
    " test/test_utils.py"
)
result = subprocess.run(cmd_str, shell=True)
# I'm not 100% sure this is necessary to pass the subprocess return code
# mostly, I want CI to properly fail if a test failed
sys.exit(result.returncode)
