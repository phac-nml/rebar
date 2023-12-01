# Run

`rebar` begins by comparing a query sequence to the dataset populations in order to find its best match. The best match is simplify defined as the population with the greatest number of shared mutations, and the least number of conflicting bases. Sites with missing data ("N") and deletions ("-") are ignored in this calculation. The best match represents the primary parent of the query sequence.

If a sequence's best match had mutational conflicts, `rebar` will search for secondary parents (recombination) by testing four different recombination hypotheses:

1. Non-Recombinant
1. Designated Recombinant (using known parents from the dataset)
1. Recursive Recombinant (allowing parents to be recombinants themselves)
1. Non-Recursive Recombinant (not allowing parents to be recombinants)

(To be continued!)
