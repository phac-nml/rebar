# Run

> **Tip**: The inner workings of the `rebar` algorithm can be exposed by including the `--verbosity debug` flag!

`rebar` begins by comparing a query sequence to the dataset populations in order to find its best match. The best match is simplify defined as the population with the greatest number of shared mutations, and the least number of conflicting bases. Sites with missing data ("N") and deletions ("-") are ignored in this calculation. The best match represents the primary parent of the query sequence.

`rebar` then proceeds to search for...

---

You can see the inner-workings of the `rebar` algorithm by using `--verbosity debug`. Let's test this on the SARS-CoV-2 recombinant `XCC` which has known parents `XBB.1.9.1` and `CH.1.1.1`.

```bash
rebar run \
    --dataset-dir dataset/sars-cov-2/2023-11-30  \
    --populations "XCC" \
    --output-dir output/example/debug \
    --verbosity debug
```

The debugging output will report detailed information on dataset searches for the primary parent (best match/conensus population). In addition, it will search for secondary parents (recombination) by testing four different recombination hypotheses:

1. Non-Recombinant
1. Designated Recombinant (using known parents from the dataset)
1. Recursive Recombinant (allowing parents to be recombinants themselves)
1. Non-Recursive Recombinant (not allowing parents to be recombinants)

The best match/primary parent is found to be for `XCC` is... itself, excellent! `XCC` is a known recombinant, so that rules out **Hypothesis \#1**.

Since `XCC` is a known recombinant, `rebar` will evaluate **Hypothesis 2**. The primary parent search will be redone focusing exclusively on designated parents (`XBB.1.9.1` and `CH.1.1.1`).

search for a secondary parent in designated parents for **Hypothesis 2** (`XBB.1.9.1` and `CH.1.1.1`.)

Since we're using a dataset population, `XCC` is an exact match to itself. So there are no mutational conflicts that need to be explained by recombination/secondary parents. No evidence for **Hypothesis #3** is found.

`rebar` will then search for a secondary parent among all possible populations. However since `XCC` is an exact match to itself, there are no mutational conflicts to

Since `XCC` is a known recombinant, `rebar` will search for a secondary parent in designated parents for **Hypothesis 2** (`XBB.1.9.1` and `CH.1.1.1`.)
`rebar` will only find evidence to support Hypotheses \#2 (designated) and \#4 (recursive).
