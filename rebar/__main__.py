#!/usr/bin/env python3

import os
from datetime import datetime
from multiprocess import cpu_count
import sys
from .utils import create_logger
from rebar import make_parser

"""
Stub function and module used as a setuptools entry point.
Based on augur and treetime's __main__.py and setup.py
"""


def main():
    parser = make_parser()
    params = parser.parse_args()

    # If params doesn't have a log object, this is the help subcommand
    if not hasattr(params, "log"):
        return_code = params.func(params)

    # Otherwise, run an actual analysis subcommand
    else:
        # Create log
        params.logger = create_logger(params.log)

        # Create output directory if it doesn't exist
        if hasattr(params, "output"):
            params.outdir = os.path.dirname(params.output)
        if not os.path.exists(params.outdir) and params.outdir != "":
            params.logger.info(
                str(datetime.now()) + "\tCreating output directory: " + params.outdir
            )
            os.mkdir(params.outdir)

        # Initialize system resources for multiprocessing
        available_cpus = cpu_count()
        if hasattr(params, "threads"):
            if params.threads > available_cpus:
                params.threads = available_cpus
            # Only print this in subcommands that use multiprocessing
            # ex. not in subs or tree, so as to not confuse the user
            # as to whether they can or can't use multiple threads
            params.logger.info(
                str(datetime.now())
                + "\tUsing {} CPUs out of {} available.".format(
                    params.threads, available_cpus
                )
            )
        else:
            params.threads = 1

        return_code = params.func(params)
    sys.exit(return_code)


if __name__ == "__main__":
    main()
