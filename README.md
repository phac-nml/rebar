# rebar

[![All Contributors](https://img.shields.io/badge/all_contributors-11-orange.svg?style=flat-square)](#credits)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://github.com/phac-nml/rebar/blob/master/LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/phac-nml/rebar.svg)](https://github.com/phac-nml/rebar/issues)
[![Build CI](https://github.com/phac-nml/rebar/actions/workflows/build.yaml/badge.svg)](https://github.com/phac-nml/rebar/actions/workflows/build.yaml)
[![Test CI](https://github.com/phac-nml/rebar/actions/workflows/test.yaml/badge.svg)](https://github.com/phac-nml/rebar/actions/workflows/test.yaml)
[![Latest CI](https://github.com/phac-nml/rebar/actions/workflows/latest.yaml/badge.svg)](https://github.com/phac-nml/rebar/actions/workflows/latest.yaml)

**RE**combination **BAR**code detector

## Why rebar?

`rebar` is a command-line application that _detects_ and _visualizes_ recombination between genomic sequences. It follows the [PHA4GE Guidance for Detecting and Characterizing SARS-CoV-2 Recombinants](https://github.com/pha4ge/pipeline-resources/blob/main/docs/sc2-recombinants.md) which outlines three steps:

1. Assess the genomic evidence for recombination.
1. Identify the breakpoint coordinates and parental regions.
1. Classify sequences as _designated_ or _novel_ recombinant lineages.

![A plot of the breakpoints and parental regions for the recombinant SARS-CoV-2 lineage XBB.1.16. At the top are rectangles arranged side-by-side horizontally. These are colored and labelled by each parent (ex. BJ.1., CJ.1) and are intepreted as reading left to right, 5' to 3'. Below these regions are genomic annotations, which show the coordinates for each gene. At the bottom are horizontal tracks, where each row is a sample, and each column is a mutation. Mutations are colored according to which parent the recombination region derives from.](assets/images/XBB_BJ.1_CJ.1_22897-22941.png)

## Install

```bash
wget -O rebar https://github.com/phac-nml/rebar/releases/download/v0.1.3/rebar-x86_64-unknown-linux-musl
./rebar --help
```

- Please see the [install](docs/install.md) docs for Windows, macOS, Docker, Singularity, and Conda.
- Please see the [compile](docs/compile.md) docs for those interested in source compilation.

## Usage

1. Preview pre-built datasets.

    ```bash
    rebar dataset list
    ```

1. Download a pre-built dataset, version-controlled to a specific date (try any date!).

    ```bash
    rebar dataset download --name sars-cov-2 --tag 2023-11-17 --output-dir dataset/sars-cov-2/2023-11-17
    ```

1. Detect recombination in dataset populations.

    ```bash
    rebar run \
      --dataset-dir dataset/sars-cov-2/2023-11-17  \
      --populations "AY.4.2*,BA.5.2,XBC.1.6*,XBB.1.5.1,XBL" \
      --output-dir example1
    ```

1. Plot breakpoints and parental regions.

    ```bash
    rebar plot --dataset-dir dataset/sars-cov-2/2023-11-17 --output-dir example1
    ```

Please see the [examples](docs/examples.md) docs for more tutorials including:

- Using an alignment of genomes as input.
- Performing a 'knockout' experiment.
- Validating all populations in a dataset.
- Running a custom dataset.

## Output

### Plots

Visualization of substitutions, parental origins, and breakpoints: `<output-dir>/plots/`. Please see the [Why rebar?](#why-rebar) section for an example!

### Table

A linelist summary of results: `<output-dir>/linelist.tsv`

|strain               |validate|validate_details|population|recombinant|parents  |breakpoints|edge_case|unique_key               |regions                          |substitutions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |genome_length|dataset_name|dataset_tag|cli_version|
|:--------------------|:-------|:---------------|:---------|:----------|:--------|:----------|:--------|:------------------------|:--------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------|:-----------|:----------|:----------|
|population_XBB.1.16  |pass    |                |XBB.1.16  |XBB        |BJ.1,CJ.1|22897-22941|false    |XBB_BJ.1_CJ.1_22897-22941|261-22896\|BJ.1,22942-29118\|CJ.1|A405G,T670G,C2790T,C3037T,G4184A,C4321T,C9344T,A9424G,C9534T,C9866T,C10029T,C10198T,G10447A,C10449A,C12880T,C14408T,G15451A,C15714T,C15738T,T15939C,T16342C,C17410T,T17859C,A18163G,C19955T,A20055G,C21618T,T21810C,G21987A,C22000A,C22109G,T22200A,G22577C,G22578A,G22599C,C22664A,C22674T,T22679C,C22686T,A22688G,G22775A,A22786C,G22813T,T22882G,G22895C,T22896C\|BJ.1;T22942G,T23018C,T23019C,T23031C,C25416T,A26275G\|CJ.1;C11750T,C11956T,T12730A,A14856G,G18703T,A19326G,A22101T,G22317T,C22995G,A22995C,G27915T,T28297C,A28447G\|private        |29903        |sars-cov-2  |2023-11-17 |0.1.0      |
|population_XBB.1.16.1|pass    |                |XBB.1.16.1|XBB        |BJ.1,CJ.1|22897-22941|false    |XBB_BJ.1_CJ.1_22897-22941|261-22896\|BJ.1,22942-29118\|CJ.1|A405G,T670G,C2790T,C3037T,G4184A,C4321T,C9344T,A9424G,C9534T,C9866T,C10029T,C10198T,G10447A,C10449A,C12880T,C14408T,G15451A,C15714T,C15738T,T15939C,T16342C,C17410T,T17859C,A18163G,C19955T,A20055G,C21618T,T21810C,G21987A,C22000A,C22109G,T22200A,G22577C,G22578A,G22599C,C22664A,C22674T,T22679C,C22686T,A22688G,G22775A,A22786C,G22813T,T22882G,G22895C,T22896C\|BJ.1;T22942G,T23018C,T23019C,T23031C,C25416T,A26275G\|CJ.1;C11750T,C11956T,T12730A,A14856G,G18703T,A19326G,A22101T,G22317T,C22995G,A22995C,C23202T,G27915T,T28297C,A28447G\|private|29903        |sars-cov-2  |2023-11-17 |0.1.0      |

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

Thanks go to the following people, who participated in the development of `rebar` and [ncov-recombinant](https://github.com/ktmeaton/ncov-recombinant):

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
    <td align="center"><a href="https://www.inspq.qc.ca/en/auteurs/2629/all"><img src="https://i1.rgstatic.net/ii/profile.image/278724097396748-1443464411327_Q128/Benjamin-Delisle.jpg?s=100" width="100px;" alt=""/><br /><sub><b>Benjamin Delisle</b></sub></a><br /><a href="https://github.com/phac-nml/rebar/issues?q=author%3Abenjamindeslisle" title="Bug reports">🐛</a> <a href="https://github.com/phac-nml/rebar/commits?author=benjamindeslisle" title="Tests">⚠️</a></td>  
    <td align="center"><a href="https://ca.linkedin.com/in/dr-vani-priyadarsini-ikkurti-4a2ab676"><img src="https://media-exp1.licdn.com/dms/image/C5603AQHaG8Xx4QLXSQ/profile-displayphoto-shrink_200_200/0/1569339145568?e=2147483647&v=beta&t=3WrvCciW-x8J3Aw4JHGrWOpuqiikrrGV2KsDaISnHIw" width="100px;" alt=""/><br /><sub><b>Vani Priyadarsini Ikkurthi</b></sub></a><br /><a href="https://github.com/phac-nml/rebar/issues?q=author%3Avanipriyadarsiniikkurthi" title="Bug reports">🐛</a> <a href="https://github.com/phac-nml/rebar/commits?author=vanipriyadarsiniikkurthi" title="Tests">⚠️</a></td>
    <td align="center"><a href="https://ca.linkedin.com/in/mark-horsman-52a14740"><img src="https://ui-avatars.com/api/?name=Mark+Horsman?s=100" width="100px;" alt=""/><br /><sub><b>Mark Horsman</b></sub></a><br /><a href="#ideas-markhorsman" title="Ideas, Planning, & Feedback">🤔</a> <a href="#design-markhorsman" title="Design">🎨</a></td>
    <td align="center"><a href="https://github.com/dfornika"><img src="https://avatars.githubusercontent.com/u/145659?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dan Fornika</b></sub></a><br /><a href="#ideas-dfornika" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/phac-nml/rebar/commits?author=dfornika" title="Tests">⚠️</a></td>
    <td align="center"><img src="https://ui-avatars.com/api/?name=Tara+Newman?s=100" width="100px;" alt=""/><br /><sub><b>Tara Newman</b></sub><br /><a href="#ideas-TaraNewman" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/phac-nml/rebar/commits?author=TaraNewman" title="Tests">⚠️</a></td>  
  </tr>
    <td align="center">
      <a href="https://github.com/TheZetner"><img src="https://avatars.githubusercontent.com/u/11616351?v=4s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Adrian Zetner</b></sub>
      </a>
      <br />
      <a href="" title="Code Review">🔣</a>
      <a href="r" title="Ideas">🤔</a>
    </td>
    <td align="center">
      <a href="https://github.com/ConnorChato"><img src="https://avatars.githubusercontent.com/u/24962136?v=4?s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Connor Chato</b></sub>
      </a>
      <br />
      <a href="" title="Code Review">🔣</a>
      <a href="r" title="Ideas">🤔</a>
    </td>
    <td align="center">
      <a href="https://github.com/mattheww95"><img src="https://avatars.githubusercontent.com/u/76452933?v=4?s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Matthew Wells</b></sub>
      </a>
      <br />
      <a href="r" title="Cross-Platoform Compilation">📦</a>
    </td>  
    <td align="center">
      <a href="https://github.com/AndreaTy"><img src="https://ui-avatars.com/api/?name=AndreaTyler?s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Andrea Tyler</b></sub>
      </a>
      <br />
      <a href="" title="Code Review">🔣</a>
    </td>
  <tr>

  </tr>
</table>
