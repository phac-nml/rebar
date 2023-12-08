# Development

- Add install documentation for conda.

## CLI

- Parse `--populations` and `--alignment` as a stream, writing to `linelist.tsv` in realtime.
- Add progress bars for parsing the `--alignment`.
- TBD: Add progress bars for exporting barcodes.

## Library

- Add parameter `recombination` to `phylogeny::get_ancestors` and `phylogeny::get_descendants`.
  - Controls whether we want to count recombination events.
- Remove function `phylogeny::get_names`.
  - Unnecessary and unused, now that I figured out how `petgraph::visit::IntoNodeReferences` works.
- Remove unnecessary parameters from struct `Recombination`: `sequence` and `genome_length`.
  - This also makes `Recombination` lifetime-free, with no more dependency on `seq`.
- Improve function `phylogeny::get_common_ancestor`.
  - This was a significant bottleneck before.
- Update function `dataset::expand_populations`.
  - Differentiate between wildcard `*` (all descendants) and `*-r` (all descendants excluding recombination).
- Switch from external crate `bio` to `noodles`.
  - `noodles` cuts down on dependencies by removing 70 external crates!
  - `noodles` will also help us prepare us for Issue #14, as it has VCF.
- `export::linelist` add new parameter `Sequence` to `results`:
  - results: &Vec<(Sequence, SearchResult, Recombination)>
