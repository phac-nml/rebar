## Compilation

For all operating systems, you will need to:

- Install the Rust compiler: [official](https://doc.rust-lang.org/cargo/getting-started/installation.html) or [conda](https://anaconda.org/conda-forge/rust).
- Clone the GitHub repository:

    ```bash
    git clone https://github.com/phac-nml/rebar.git
    cd rebar
    ```

For Linux and Windows, install [Docker](https://docs.docker.com/engine/install/).

1. Clone the [rebar](https://github.com/phac-nml/rebar) repository.

    ```bash
    git clone https://github.com/phac-nml/rebar.git
    cd rebar
    ```

### Linux

1. Install [Docker](https://docs.docker.com/engine/install/).

1. Pull the [Rust compiler](https://rustup.rs/).

1. Install the [Rust compiler](https://rustup.rs/).

1. Clone the [rebar](https://github.com/phac-nml/rebar) repository.

    ```bash
    git clone https://github.com/phac-nml/rebar.git
    cd rebar
    ```

1. Install [cross](https://github.com/cross-rs/cross) to compile from Docker images.

    ```bash
    cargo install cross
    ```

1. Compile.

    ```bash
    cross build --release --all-features --target x86_64-unknown-linux-musl
    ```

    - We use `musl` rather than `gnu` to generate standalone binaries that are not dependent on system C libraries (GLIBC). This improves compatibility with "older" systems (ex. CentOS 7). Although we could use the tag `main-centos` if we really want to... (ex. image: `ghcr.io/cross-rs/x86_64-unknown-linux-gnu:main-centos`).`

1. Test.

    ```bash
    target/x86_64-unknown-linux-musl/release/rebar
    ```

### Windows

1. Install the [Rust compiler](https://rustup.rs/).

    - Manual
    - Continue without C++ libraries: Y
    - Customize installation
    - Default host triple: x86_64-pc-windows-gnu
    - Default toolchain: stable
    - Profile: minimal
    - Modify Path: Y
    - Proceed with intallation: Y

    > **Note**: If the screen output hangs for a long time on `info: installing component ...`, press `<ENTER>` a couple times.

1. Install [Docker](https://docs.docker.com/engine/install/) and start the Docker engine.

1. Install [cross](https://github.com/cross-rs/cross).

    ```bash
    cargo install cross
    ```

1. Compile.

    ```bash
    cross build --release --all-features --target x86_64-pc-windows-gnu
    ```
