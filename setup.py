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

with open("requirements.txt", "r") as r:
    require_list = r.read().strip().split("\n")

setup(
    name="bio-rebar",
    version=get_version(),
    author="Katherine Eaton",
    author_email="katherine.eaton@phac-aspc.gc.ca",
    description=("REcombination BARcode detection for SARS-CoV-2."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="SARS-CoV-2, recombination",
    url="https://github.com/phac-nml/rebar",
    packages=["rebar"],
    install_requires=require_list,
    entry_points={
        "console_scripts": [
            "rebar = rebar.__main__:main",
        ]
    },
)
