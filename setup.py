from setuptools import setup


def get_version():
    version = "0.0.0"
    with open("rebar/__init__.py") as ifile:
        for line in ifile:
            if line[:7] == "version":
                version = line.split("=")[-1].strip()[1:-1]
                break
    return version


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="ncov-rebar",
    version=get_version(),
    author="Katherine Eaton",
    author_email="katherine.eaton@phac-aspc.gc.ca",
    description=("SARS-CoV-2 recombination detection from lineage barcodes."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="SARS-CoV-2, recombination",
    url="https://github.com/phac-nml/rebar",
    packages=["rebar"],
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.4.1",
        "pango_aliasor>=0.2.2",
        "snipit>=1.0.7",
    ],
    entry_points={
        "console_scripts": [
            "rebar = rebar.__main__:main",
        ]
    },
)
