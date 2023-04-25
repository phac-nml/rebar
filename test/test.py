#!/usr/bin/env python3

import subprocess
import sys

cmd_str = (
    "python -m coverage run -m pytest --cov=rebar --cov-report=html --cov-report=xml"
    " test/test_utils.py"
    " test/test_wrappers.py"
)
result = subprocess.run(cmd_str, shell=True)
# I'm not 100% sure this is necessary to pass the subprocess return code
# mostly, I want CI to properly fail if a test failed
sys.exit(result.returncode)
