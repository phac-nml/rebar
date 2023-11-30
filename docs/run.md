# Run

> **Tip**: The inner workings of the `rebar` algorithm can be exposed by including the `--verbosity debug` flag!

`rebar` begins by comparing a query sequence to the dataset populations in order to find its best match. The best match is simplify defined as the population with the greatest number of shared mutations, and the least number of conflicting bases. Sites with missing data ("N") and deletions ("-") are ignored in this calculation. The best match represents the primary parent of the query sequence.

`rebar` then proceeds to search for
