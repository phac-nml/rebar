# rebar

[![All Contributors](https://img.shields.io/badge/all_contributors-11-orange.svg?style=flat-square)](#credits)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://github.com/phac-nml/rebar/blob/master/LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/phac-nml/rebar.svg)](https://github.com/phac-nml/rebar/issues)
[![Build CI](https://github.com/phac-nml/rebar/actions/workflows/build.yaml/badge.svg)](https://github.com/phac-nml/rebar/actions/workflows/build.yaml)
[![Test CI](https://github.com/phac-nml/rebar/actions/workflows/test.yaml/badge.svg)](https://github.com/phac-nml/rebar/actions/workflows/test.yaml)
[![Latest CI](https://github.com/phac-nml/rebar/actions/workflows/latest.yaml/badge.svg)](https://github.com/phac-nml/rebar/actions/workflows/latest.yaml)

`rebar` is a command-line application that _detects_ and _visualizes_ recombination between genomic sequences. It follows the [PHA4GE Guidance for Detecting and Characterizing SARS-CoV-2 Recombinants](https://github.com/pha4ge/pipeline-resources/blob/main/docs/sc2-recombinants.md) which outlines three steps:

1. Assess the genomic evidence for recombination.
1. Identify the breakpoint coordinates and parental regions.
1. Classify sequences as _designated_ or _novel_ recombinant lineages.

![plot_XBB.1.16](images/XBB_BJ.1_CJ.1_22897-22941.png)

## Install

Download the [latest release](https://github.com/phac-nml/rebar/releases/download/v0.1.0/rebar-x86_64-unknown-linux-gnu).

  ```bash
  wget -O rebar https://github.com/phac-nml/rebar/releases/download/v0.1.0/rebar-x86_64-unknown-linux-gnu

  ./rebar --help
  ```

### Conda \*\*Coming Soon!\*\*

### Source

1. Clone repository: `git clone https://github.com/phac-nml/rebar.git && cd rebar`
1. Install the rust compiler: [official](https://doc.rust-lang.org/cargo/getting-started/installation.html) or [conda](https://anaconda.org/conda-forge/rust).
1. Compile: `cargo build --release --all-features --target x86_64-unknown-linux-gnu`
    - Output binary: `target/x86_64-unknown-linux-gnu/release/rebar`

## Usage

Download a sars-cov-2 dataset, version-controlled to a specific date.

  ```bash
  rebar dataset download \
    --name sars-cov-2 \
    --tag 2023-11-17 \
    --output-dir dataset/sars-cov-2/2023-11-17
  ```

- `--tag` can be any date (YYYY-MM-DD)!

### Example 1

1. Detect recombination in user-specified populations.

    ```bash
    rebar run \
      --dataset-dir dataset/sars-cov-2/2023-11-17  \
      --populations "AY.4.2*,BA.5.2,XBC.1.6*,XBB.1.5.1,XBL" \
      --output-dir example1
    ```

    - `--populations` can include any sequence name found in the dataset `populations.fasta`. For sars-cov-2, sequence names are the designated lineages.
    - The wildcard character ("\*") will include the lineage and all its descendants.
    - **NOTE**: If using "\*", make sure to use quotes (ex. `--lineages "XBC*,XBB.1.16*"`)!

1. Plot breakpoints and parental regions.

    ```bash
    rebar plot --dataset-dir dataset/sars-cov-2/2023-11-17 --output-dir example1
    ```

### Example 2

1. Detect recombination from an input alignment.

    ```bash
    rebar run \
      --dataset dataset/sars-cov-2/2023-11-17 \
      --alignment data/example2.fasta \
      --output-dir example2
    ```

    - The `--alignment` must be aligned to the same reference as in the dataset `reference.fasta` (we strongly recommend [nextclade](https://clades.nextstrain.org/)).

1. Plot breakpoints and parental regions.

    ```bash
    rebar plot --dataset-dir dataset/sars-cov-2/2023-11-17 --output-dir example2
    ```

## Validate

Run `rebar` on all populations in the dataset, and validate against the expected results.

```bash
rebar run \
  --dataset-dir dataset/sars-cov-2/latest \
  --output-dir validate \
  --populations "*" \
  --threads 4
```

## Output

### Plots

Visualization of substitutions, parental origins, and breakpoints.

### Table

A linelist summary of detection results.

|strain               |validate|validate_details|population|recombinant|parents  |breakpoints|edge_case|unique_key               |regions                          |private|diagnostic|genome_length|dataset_name|dataset_tag|cli_version|
|:--------------------|:-------|:---------------|:---------|:----------|:--------|:----------|:--------|:------------------------|:--------------------------------|:------|:---------|:------------|:-----------|:----------|:----------|
|XBB.1.16  |pass    |                |XBB.1.16  |XBB        |BJ.1,CJ.1|22897-22941|false    |XBB_BJ.1_CJ.1_22897-22941|261-22896\|BJ.1,22942-29118\|CJ.1|       |NA        |29903        |sars-cov-2  |2023-11-17 |0.1.0      |
|XBB.1.16.1|pass    |                |XBB.1.16.1|XBB        |BJ.1,CJ.1|22897-22941|false    |XBB_BJ.1_CJ.1_22897-22941|261-22896\|BJ.1,22942-29118\|CJ.1|       |NA        |29903        |sars-cov-2  |2023-11-17 |0.1.0      |

## Credits

[rebar](https://github.com/phac-nml/rebar) is built and maintained by [Katherine Eaton](https://ktmeaton.github.io/) at the [National Microbiology Laboratory (NML)](https://github.com/phac-nml) of the Public Health Agency of Canada (PHAC).

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)). Contributions of any kind welcome!

<table>
  <tr>
    <td align="center"><a href="https://ktmeaton.github.io"><img src="https://s.gravatar.com/avatar/0b9dc28b3e64b59f5ce01e809d214a4e?s=80" width="100px;" alt=""/><br /><sub><b>Katherine Eaton</b></sub></a><br /><a href="https://github.com/phac-nml/rebar/commits?author=ktmeaton" title="Code">💻</a> <a href="https://github.com/phac-nml/rebar/commits?author=ktmeaton" title="Documentation">📖</a> <a href="#design-ktmeaton" title="Design">🎨</a> <a href="#ideas-ktmeaton" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-ktmeaton" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-ktmeaton" title="Maintenance">🚧</a></td>
  </tr>
</table>

Special thanks go to the following people, who are instrumental to the design and data sources in `rebar`:

- Lena Schimmel ([@lenaschimmel](https://github.com/lenaschimmel)) for the original concept of a barcode-scanning recombinant detector with [sc2rf](https://github.com/lenaschimmel/sc2rf).
- Cornelius Roemer ([@corneliusroemer](https://github.com/corneliusroemer)) for the [designated lineages](https://github.com/cov-lineages/pango-designation), [consensus sequences](https://github.com/yatisht/usher), and [Nextclade barcodes](https://raw.githubusercontent.com/corneliusroemer/pango-sequences/main/data/pango-consensus-sequences_summary.json).
- Josh Levy ([@joshuailevy](https://github.com/andersen-lab/Freyja-data)) and the [Andersen Lab](https://github.com/andersen-lab) for the [UShER barcodes](https://github.com/yatisht/usher) from [Freyja](https://github.com/andersen-lab/Freyja).
- Richard Neher ([@rneher](https://github.com/rneher)) and the [Neher Lab](https://github.com/neherlab) for python package structure, specifically [treetime](https://github.com/neherlab/treetime).

<table>
  <tr>
    <td align="center">
      <a href="https://github.com/lenaschimmel"><img src="https://avatars.githubusercontent.com/u/1325019?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Lena Schimmel</b></sub>
      </a>
      <br />
      <a href="https://github.com/lenaschimmel/sc2rf" title="Ideas: sc2rf">🤔</a>
    </td>
    <td align="center">
      <a href="https://github.com/corneliusroemer">
        <img src="https://avatars.githubusercontent.com/u/25161793?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Cornelius Roemer</b></sub>
      </a>
      <br />
      <a href="https://github.com/cov-lineages/pango-designation" title="Data: Lineage Designations">🔣</a>
      <a href="https://github.com/corneliusroemer/pango-sequences" title="Data: Consensus Sequences">🔣</a>
      <a href="https://github.com/corneliusroemer/pango-sequences" title="Data: Nextclade Barcodes">🔣</a>
    </td>
    <td align="center">
      <a href="https://github.com/joshuailevy">
      <img src="https://avatars.githubusercontent.com/u/19437463?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Josh Levy</b></sub>
      </a>
      <br />
      <a href="https://github.com/andersen-lab/Freyja-data" title="Data: UShER Barcodes">🔣</a>
    </td>
    <td align="center">
      <a href="https://github.com/rneher">
      <img src="https://avatars.githubusercontent.com/u/8379168?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Richard Neher</b></sub>
      </a>
      <br />
      <a href="https://github.com/neherlab/treetime" title="Ideas: Treetime">🤔</a>
    </td>  
  </tr>
</table>

Thanks go to the following people, who participated in the development of [ncov-recombinant](https://github.com/ktmeaton/ncov-recombinant), which `rebar` is based on:

<table>
  <tr>
    <td align="center">
      <a href="https://github.com/yatisht"><img src="https://avatars.githubusercontent.com/u/34664884?v=4s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Yatish Turakhia</b></sub>
      </a>
      <br />
      <a href="https://github.com/yatisht/usher" title="Data: UShER">🔣</a>
      <a href="https://github.com/yatisht/usher" title="Ideas: UShER">🤔</a>
    </td>
    <td align="center">
      <a href="https://github.com/AngieHinrichs"><img src="https://avatars.githubusercontent.com/u/186983?v=4?v=4s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Angie Hinrichs</b></sub>
      </a>
      <br />
      <a href="https://github.com/yatisht/usher" title="Data: UShER">🔣</a>
      <a href="https://github.com/yatisht/usher" title="Ideas: UShER">🤔</a>
    </td>
    <td align="center"><a href="https://www.inspq.qc.ca/en/auteurs/2629/all"><img src="https://i1.rgstatic.net/ii/profile.image/278724097396748-1443464411327_Q128/Benjamin-Delisle.jpg?s=100" width="100px;" alt=""/><br /><sub><b>Benjamin Delisle</b></sub></a><br /><a href="https://github.com/phac-nml/rebar/issues?q=author%3Abenjamindeslisle" title="Bug eports">🐛</a> <a href="https://github.com/phac-nml/rebar/commits?author=benjamindeslisle" title="Tests">⚠️</a></td>  
    <td align="center"><a href="https://ca.linkedin.com/in/dr-vani-priyadarsini-ikkurti-4a2ab676"><img src="https://media-exp1.licdn.com/dms/image/C5603AQHaG8Xx4QLXSQ/profile-displayphoto-shrink_200_200/0/1569339145568?e=2147483647&v=beta&t=3WrvCciW-x8J3Aw4JHGrWOpuqiikrrGV2KsDaISnHIw" width="100px;" alt=""/><br /><sub><b>Vani Priyadarsini Ikkurthi</b></sub></a><br /><a href="https://github.com/phac-nml/rebar/issues?q=author%3Avanipriyadarsiniikkurthi" title="Bug reports">🐛</a> <a href="https://github.com/phac-nml/rebar/commits?author=vanipriyadarsiniikkurthi" title="Tests">⚠️</a></td>
    <td align="center"><a href="https://ca.linkedin.com/in/mark-horsman-52a14740"><img src="https://ui-avatars.com/api/?name=Mark+Horsman?s=100" width="100px;" alt=""/><br /><sub><b>Mark Horsman</b></sub></a><br /><a href="#ideas-markhorsman" title="Ideas, Planning, & Feedback">🤔</a> <a href="#design-markhorsman" title="Design">🎨</a></td>
    <td align="center"><a href="https://github.com/dfornika"><img src="https://avatars.githubusercontent.com/u/145659?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dan Fornika</b></sub></a><br /><a href="#ideas-dfornika" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/phac-nml/rebar/commits?author=dfornika" title="Tests">⚠️</a></td>
    <td align="center"><img src="https://ui-avatars.com/api/?name=Tara+Newman?s=100" width="100px;" alt=""/><br /><sub><b>Tara Newman</b></sub><br /><a href="#ideas-TaraNewman" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/phac-nml/rebar/commits?author=TaraNewman" title="Tests">⚠️</a></td>  
  </tr>  

</table>
