# This code is executed when imported as a module

# -----------------------------------------------------------------------------
# Version

version = "0.1.0"

# -----------------------------------------------------------------------------
# Errors


class RebarError(Exception):
    """
    RebarError class
    Parent class for more specific errors
    Raised when rebar is used incorrectly in a known way.
    """

    pass


class RebarUnknownError(Exception):
    """
    RebarUnknownError class
    Raised when rebar is used incorrectly in a unknown way.
    """

    pass


class MissingDataError(RebarError):
    """MissingDataError class raised when X is missing."""

    pass


# -----------------------------------------------------------------------------
# Module functions

from .argument_parser import make_parser
