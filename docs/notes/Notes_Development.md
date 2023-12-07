# Development

- Add install documentation for conda.

## CLI

- Parse `--populations` and `--alignment` as a stream, writing to `linelist.tsv` in realtime.
- Add progress bars for parsing the `--alignment`.
- TBD: Add progress bars for exporting barcodes.

## Library

- Remove unnecessary parameters from struct `Recombination`: `sequence` and `genome_length`.
  - This also makes `Recombination` lifetime-free, with no more dependency on `seq`.
- Improve efficiency of `Phylogeny` method `get_common_ancestor`.
  - This was a significant bottleneck before.
- Switch from external crate `bio` to `noodles`.
  - `noodles` cuts down on dependencies by removing 70 external crates!
  - `noodles` will also help us prepare us for Issue #14, as it has VCF.
