# Install

## Direct

### Linux

```bash
wget -O rebar https://github.com/phac-nml/rebar/releases/latest/download/rebar-x86_64-unknown-linux-musl
./rebar --help
```

### macOS

```bash
wget -O rebar https://github.com/phac-nml/rebar/releases/latest/download/rebar-x86_64-apple-darwin
./rebar --help
```

#### Windows

```powershell
wget "https://github.com/phac-nml/rebar/releases/latest/download/rebar-x86_64-pc-windows-gnu.exe" -OutFile rebar.exe
.\rebar.exe --help
```

## Conda

> **Coming Soon!** The conda install option will be available when this page is live: https://anaconda.org/bioconda/rebar

```bash
conda create -n rebar rebar
conda activate rebar
rebar --help
```

## Docker

### GitHub Release

1. Download and load the Docker image from the releases page.

    ```bash
    wget -O rebar-docker.tar https://github.com/phac-nml/rebar/releases/latest/download/rebar-docker.tar
    docker load --input rebar-docker.tar
    ```

1. Run the image.

    ```bash
    docker run -v .:/rebar ghcr.io/phac-nml/rebar:latest \
    rebar --help
    ```

    > *TIP*: Use `-v .:/rebar` to connect your local directory (`.`) to the docker image. That way when you run `rebar dataset download`, it will download permanently onto your system, not temporarily inside the image.

### GitHub Container Registry

At the moment (2023-12-04), the Github Container Registry is locked to 'private'. Only `phac-nml` members can access. We are working on this!

1. Authenticate with the GitHub Container Registry.

    ```bash
    export GITHUB_USERNAME='<Github Username>'
    export GITHUB_TOKEN='<Github Personal Access Token>'
    echo $GITHUB_TOKEN | docker login ghcr.io -u $GITHUB_USERNAME --password-stdin
    ```

1. Run the image.

    ```bash
    docker run -v .:/rebar ghcr.io/phac-nml/rebar:latest \
    rebar --help
    ```

## Singularity

### GitHub Release

1. Download and build a singularity image from the docker image.

    ```bash
    wget -O rebar-docker.tar https://github.com/phac-nml/rebar/releases/latest/download/rebar-docker.tar
    singularity build rebar_latest.sif docker-archive://rebar-docker.tar
    ```

1. Run the image

    ```bash
    singularity run --home $(pwd) rebar_latest.sif \
    rebar --help
    ```

    > **Tip**: Use `--home $(pwd)` to have `rebar` run in your current working directory.

### Container Registry

At the moment (2023-12-04), the Github Container Registry is locked to 'private'. Only `phac-nml` members can access. We are working on this!

1. Download the image from the GitHub Container Registry.

    ```bash
    export SINGULARITY_DOCKER_USERNAME='<Github Username>'
    export SINGULARITY_DOCKER_PASSWORD='<Github Personal Access Token>'
    singularity pull docker://ghcr.io/phac-nml/rebar:latest
    ```

1. Run the image

    ```bash
    singularity run --home $(pwd) rebar_latest.sif \
    rebar --help
    ```
