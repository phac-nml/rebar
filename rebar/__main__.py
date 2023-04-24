#!/usr/bin/env python3

import os
from datetime import datetime
from multiprocess import cpu_count
import sys
from .utils import create_logger
from rebar import make_parser, RebarError
import timeit

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
        sys.exit(return_code)

    # Check for at least one output type specified
    if hasattr(params, "output_all"):
        if (
            not params.output_all
            and not params.output_fasta
            and not params.output_barcode
            and not params.output_plot
            and not params.output_tsv
            and not params.output_yaml
        ):
            raise SystemExit(
                RebarError(
                    "RebarError: At least one output type must be specified"
                    " with --output-TYPE."
                )
            )

    # Check for conflicting use of --debug and --threads
    if params.debug and params.threads > 1:
        raise SystemExit(
            RebarError(
                "RebarError: Debugging mode (--debug) and multithreading mode"
                " (--threads) are incompatible. Please specify only one or the other."
            )
        )
    
    # Check for validate mode and missing tsv output  
    if hasattr(params, "validate"):
        if params.validate and not (params.output_all or params.output_tsv):
            raise SystemExit(
                RebarError(
                    "RebarError: --validate requires either --output-all or --output-tsv."
                )
            )    

    # Otherwise, run an actual analysis subcommand
    # Create log directory if it doesn't increase
    if params.log:
        logdir = os.path.dirname(params.log)
        if not os.path.exists(logdir) and logdir != "":
            os.makedirs(logdir)

    # Create log
    params.logger = create_logger(params.log)

    # Create output directory if it doesn't exist
    if hasattr(params, "output"):
        params.outdir = os.path.dirname(params.output)
    if not os.path.exists(params.outdir) and params.outdir != "":
        params.logger.info(
            str(datetime.now()) + "\tCreating output directory: " + params.outdir
        )
        os.makedirs(params.outdir)

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

    # Reverse the no_edge_cases parameter
    if hasattr(params, "no_edge_cases"):
        params.edge_cases = True
        if params.no_edge_cases:
            params.edge_cases = False

    return_code = params.func(params)

    sys.exit(return_code)


if __name__ == "__main__":
    main()
