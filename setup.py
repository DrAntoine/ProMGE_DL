#!/usr/bin/env python
# -*-coding: utf-8 -*-

from promgedl.version import __version__

from setuptools import setup, find_packages

url = "https://github.com/DrAntoine/ProMGE_DL"
desc = "A tool to easily use the ProMGE database and associated sequences"

setup(
    name="promgedl",
    version=__version__,

    description=desc,
    url=url,
    author="Antoine Druart",
    author_email="antoine.druart@me.com",

    license="MIT",
    packages=find_packages(),

    install_requires = [
        "biopython==1.81",
        "tqdm"
    ],

    entry_points={
        "console_scripts": ["promgedl = promgedl.app:main"],
    }
)