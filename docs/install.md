# Install

## Direct

### Linux

```bash
wget -O rebar https://github.com/phac-nml/rebar/releases/download/v0.1.3/rebar-x86_64-unknown-linux-musl
./rebar --help
```

### macOS

```bash
wget -O rebar https://github.com/phac-nml/rebar/releases/download/v0.1.3/rebar-x86_64-apple-darwin
./rebar --help
```

#### Windows

```powershell
wget "https://github.com/phac-nml/rebar/releases/download/v0.1.3/rebar-x86_64-pc-windows-gnu" -OutFile rebar.exe
.\rebar.exe --help
```

## Docker

```bash
docker run -v .:/rebar ghcr.io/phac-nml/rebar:latest \
rebar --help
```

## Singularity

```bash
singularity run docker://ghcr.io/phac-nml/rebar:latest \
rebar --help
```

## Conda

\*\*Coming Soon!\*\*
