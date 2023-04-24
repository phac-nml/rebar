# rebar

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-10-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://github.com/phac-nml/rebar/blob/master/LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/phac-nml/rebar.svg)](https://github.com/phac-nml/rebar/issues)
[![Install CI](https://github.com/phac-nml/rebar/actions/workflows/install.yaml/badge.svg)](https://github.com/phac-nml/rebar/actions/workflows/install.yaml)
[![Pipeline CI](https://github.com/phac-nml/rebar/actions/workflows/pipeline.yaml/badge.svg)](https://github.com/phac-nml/rebar/actions/workflows/pipeline.yaml)

**RE**combination detection using **BAR**codes for SARS-CoV-2.

## Install

1. Source

    ```bash
    git clone https://github.com/phac-nml/rebar.git
    cd rebar
    pip install .
    ```

1. PyPI

    \*\*Coming Soon\*\*

1. Conda

    \*\*Coming Soon\*\*

## Usage

Begin by downloading and building the `rebar` dataset model with:

```bash
rebar dataset \
  --name sars-cov-2 \
  --tag latest \
  --outdir dataset/sars-cov-2-latest
```

The `dataset` subcommand performs the following tasks:

- Downloads the latest _designated_ lineages from [pango-designation](https://github.com/cov-lineages/pango-designation).
- Creates a nomenclature tree of designated lineages.
- Downloads lineage barcodes from [pango-sequences](https://github.com/corneliusroemer/pango-sequences) (Nextclade) and [Freyja-data](https://github.com/andersen-lab/Freyja-data) (UShER).

### Example 1

Detect recombination in user-specified lineages.

```bash
rebar run \
  --dataset dataset/sars-cov-2-latest \
  --lineages AY.4,BA.5.2,XD,XBB.1.5.1,XBL \
  --output-all \
  --outdir example1
```

The `--lineages` can include any designated lineage found in the dataset `alignment.fasta`.

### Example 2

Detect recombination from an input alignment:

```bash
rebar run \
  --dataset dataset/sars-cov-2-latest \
  --alignment test/alignment.fasta \
  --output-all \
  --outdir example2
```

The `--alignment` must be aligned to the same reference as in the dataset `reference.fasta` (we strongly recommend [nextclade](https://github.com/nextstrain/nextclade)).

### Example 3

Detect recombination in all designated lineages in the dataset model.

```bash
rebar run \
  --dataset dataset/sars-cov-2-latest \
  --alignment dataset/sars-cov-2-latest/alignment.fasta \
  --validate \
  --threads 8 \
  --output-all \
  --outdir example3
```

## Output

### Plot

Visualization with [snipit](https://github.com/aineniamh/snipit): `snipit_<parent_1>_<parent_2>.png`

![snipit_XD](images/snipit_XD.png)

### Table

A linelist summary of recombination detection results: `summary.tsv`

|strain   |lineage  |clade      |recombinant|parents_lineage|parents_clade|breakpoints            |regions                                            |
|:--------|:--------|:----------|:----------|:--------------|:------------|:----------------------|:--------------------------------------------------|
|AY.4     |AY.4     |21J        |False      |               |             |                       |                                                   |
|BA.5.2   |BA.5.2   |22B        |False      |               |             |                       |                                                   |
|XBB.1.5.1|XBB.1.5.1|23A        |XBB        |BJ.1,CJ.1      |21L,22D      |22897:22941            |261-22896\|BJ.1,22942-29118\|CJ.1                  |
|XBL      |XBL      |recombinant|XBL        |XBB.1.5,BA.2.75|22F,22D      |5184:16341             |3796-5183\|BA.2.75,16342-27915\|XBB.1.5            |
|XD       |XD       |recombinant|XD         |BA.1,AY.4      |21K,21J      |21988:22577,25470:25583|210-21987\|AY.4,22578-25469\|BA.1,25584-29402\|AY.4|

### Summary YAML

A detailed YAML summary of recombination detection results: `summary.yaml`

```yaml
XA:
  substitutions: C241T,T445C,T2019C,C3037T,C4999T,C6286T,C8090T,C9430T,A10323G,C13945T,C14408T,G20410A,G21255C,A23063T,C23208T,C23271A,A23403G,C23604A,C23709T,T24506G,G2491
  deletions: 21765-21770,21992-21994,28271
  missing: 1-200,29704-29903
  lineage:
    lineage: XA
    clade: recombinant
    top_lineages: XA
    top_lineages_subsample: XA
    outlier_lineages:
    barcode: C241T,T445C,T2019C,C3037T,C4999T,C6286T,C8090T,C9430T,A10323G,C13945T,C14408T,G20410A,G21255C,A23063T,C23208T,C23271A,A23403G,C23604A,C23709T,T24506G,G24914C,G
    support: C241T,T445C,T2019C,C3037T,C4999T,C6286T,C8090T,C9430T,A10323G,C13945T,C14408T,G20410A,G21255C,A23063T,C23208T,C23271A,A23403G,C23604A,C23709T,T24506G,G24914C,G
    missing:
    conflict_ref:
    conflict_alt:
    recombinant: XA
    recursive: None
    edge_case: 'False'
    ...  
```

## Credits

[rebar](https://github.com/phac-nml/rebar) is built and maintained by [Katherine Eaton](https://ktmeaton.github.io/) at the [National Microbiology Laboratory (NML)](https://github.com/phac-nml) of the Public Health Agency of Canada (PHAC).

<table>
  <tr>
    <td align="center"><a href="https://ktmeaton.github.io"><img src="https://s.gravatar.com/avatar/0b9dc28b3e64b59f5ce01e809d214a4e?s=80" width="100px;" alt=""/><br /><sub><b>Katherine Eaton</b></sub></a><br /><a href="https://github.com/phac-nml/rebar/commits?author=ktmeaton" title="Code">💻</a> <a href="https://github.com/phac-nml/rebar/commits?author=ktmeaton" title="Documentation">📖</a> <a href="#design-ktmeaton" title="Design">🎨</a> <a href="#ideas-ktmeaton" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-ktmeaton" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-ktmeaton" title="Maintenance">🚧</a></td>
  </tr>
</table>

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

- Note: Add Cornelius Roemer, neherlab (treetime)

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/nextstrain/nextclade"><img src="https://avatars.githubusercontent.com/u/22159334?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Nextstrain (Nextclade)</b></sub></a><br /><a href="#data-nextstrain" title="Data">🔣</a> <a href="#plugin-nextstrain" title="Plugin/utility libraries">🔌</a></td>
    <td align="center"><a href="https://github.com/lenaschimmel/sc2rf"><img src="https://avatars.githubusercontent.com/u/1325019?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Lena Schimmel (sc2rf)</b></sub></a><br /><a href="#plugin-lenaschimmel" title="Plugin/utility libraries">🔌</a></td>
    <td align="center"><a href="https://github.com/yatisht/usher"><img src="https://avatars.githubusercontent.com/u/34664884?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Yatish Turakhia (UShER)</b></sub></a><br /><a href="#data-yatisht" title="Data">🔣</a> <a href="#plugin-yatisht" title="Plugin/utility libraries">🔌</a></td>
    <td align="center"><a href="https://github.com/yatisht/usher"><img src="https://avatars.githubusercontent.com/u/186983?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Angie Hinrichs (UShER)</b></sub></a><br /><a href="#data-AngieHinrichs" title="Data">🔣</a> <a href="#plugin-AngieHinrichs" title="Plugin/utility libraries">🔌</a></td>
    <td align="center"><a href="https://www.inspq.qc.ca/en/auteurs/2629/all"><img src="https://i1.rgstatic.net/ii/profile.image/278724097396748-1443464411327_Q128/Benjamin-Delisle.jpg?s=100" width="100px;" alt=""/><br /><sub><b>Benjamin Delisle</b></sub></a><br /><a href="https://github.com/phac-nml/rebar/issues?q=author%3Abenjamindeslisle" title="Bug reports">🐛</a> <a href="https://github.com/phac-nml/rebar/commits?author=benjamindeslisle" title="Tests">⚠️</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://ca.linkedin.com/in/dr-vani-priyadarsini-ikkurti-4a2ab676"><img src="https://media-exp1.licdn.com/dms/image/C5603AQHaG8Xx4QLXSQ/profile-displayphoto-shrink_200_200/0/1569339145568?e=2147483647&v=beta&t=3WrvCciW-x8J3Aw4JHGrWOpuqiikrrGV2KsDaISnHIw" width="100px;" alt=""/><br /><sub><b>Vani Priyadarsini Ikkurthi</b></sub></a><br /><a href="https://github.com/phac-nml/rebar/issues?q=author%3Avanipriyadarsiniikkurthi" title="Bug reports">🐛</a> <a href="https://github.com/phac-nml/rebar/commits?author=vanipriyadarsiniikkurthi" title="Tests">⚠️</a></td>
    <td align="center"><a href="https://ca.linkedin.com/in/mark-horsman-52a14740"><img src="https://ui-avatars.com/api/?name=Mark+Horsman?s=100" width="100px;" alt=""/><br /><sub><b>Mark Horsman</b></sub></a><br /><a href="#ideas-markhorsman" title="Ideas, Planning, & Feedback">🤔</a> <a href="#design-markhorsman" title="Design">🎨</a></td>
    <td align="center"><a href="https://github.com/jbloomlab"><img src="https://avatars.githubusercontent.com/u/17679492?s=200&v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jesse Bloom Lab</b></sub></a><br /><a href="#data-jbloomlab" title="Data">🔣</a> <a href="#plugin-jbloomlab" title="Plugin/utility libraries">🔌</a></td>
    <td align="center"><a href="https://github.com/dfornika"><img src="https://avatars.githubusercontent.com/u/145659?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dan Fornika</b></sub></a><br /><a href="#ideas-dfornika" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/phac-nml/rebar/commits?author=dfornika" title="Tests">⚠️</a></td>
    <td align="center"><img src="https://ui-avatars.com/api/?name=Tara+Newman?s=100" width="100px;" alt=""/><br /><sub><b>Tara Newman</b></sub><br /><a href="#ideas-TaraNewman" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/phac-nml/rebar/commits?author=TaraNewman" title="Tests">⚠️</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
