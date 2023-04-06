#!/usr/bin/env python3

"""
Stub function and module used as a setuptools entry point.
Based on augur and treetime's __main__.py and setup.py
"""

import sys
from rebar import make_parser


def main():
    parser = make_parser()
    params = parser.parse_args()
    return_code = params.func(params)
    sys.exit(return_code)


if __name__ == "__main__":
    main()
