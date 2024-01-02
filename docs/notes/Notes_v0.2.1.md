# v0.2.1

This patch release documents the new conda installation method, resolves a plotting artifact and improves performance of the `phylogeny` function `get_common_ancestor`.

## New

- Issue #10, PR #34 | Add documentation for Conda installation.

## Fixes

- Issue #27, PR #37 | Fix legend overflow in plot.

## Changes

- Issue #23, PR #39 | Change `parsimony::from_sequence` param `coordinates` from Vector to Slice.
- Issue #28, PR #41 | Improve peformance of `phylogeny::get_common_ancestor`
- Issue #29, PR #38 | Reduce large artifact uploads in CI.
