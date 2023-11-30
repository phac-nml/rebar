# Dataset

A `rebar` dataset consists of two mandatory parts:

1. `reference.fasta`: The reference genome.

    ```text
    >Reference
    AAAAAAAAAAAAAAAAAAAA
    ```

1. `populations.fasta`: Known populations (ex. clades, lineages) aligned to the reference.

    ```text
    >A
    CCCCCCAACCCCCCCCCCCC
    >B
    TTTTTTTTTTTTTTTTTTAA
    >C
    AAGGGGGGGGGGGGGGGGGG
    >D
    CCCCCCAACCCTTTTTTTAA
    >E
    AAGCCCAACCCTTTTTTTAA
    ```

The following are optional components:

1. `phylogeny.json`: A phylogenetic graph which provides prior information about the evolutionary history. This is particularly useful if populations in `populations.fasta` are internal nodes or known recombinants.

    ```json
    {
      "graph": {
        "nodes": [ "root", "A", "B", "C"],
        "edge_property": "directed",
        "edges": [
          [ 0, 1, 1],
          [ 0, 2, 1],
          [ 1, 3, 1],
          [ 2, 3, 1]
        ]
      }
    }
    ```

1. `annotations.tsv`: A table of genome annotations to add to the plot.

    |gene |abbreviation|start|end|
    |:----|:-----------|:----|:--|
    |Gene1|g1          |1    |3  |
    |Gene2|g2          |7    |10 |
