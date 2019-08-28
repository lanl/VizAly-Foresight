import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vizaly-pat_USERNAME",
    version="0.0.1",
    author="Pascal Grosset, Chris Biwer, Jesus Pulido",
    author_email="exasky@lanl.com",
    description="Toolkit for running analysis at scale using slurm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages
    url="https://github.com/lanl/VizAly-Foresight/tree/master/Analysis/pat",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD-3 ",
        "Operating System :: OS Independent",
    ],
)