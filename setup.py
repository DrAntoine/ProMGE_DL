#!/usr/bin/env python
# -*-coding: utf-8 -*-

from promgedl.version import __version__

from setuptools import setup, find_packages

url = "https://github.com/DrAntoine/ProMGE_DL"
desc = "Download the ProMGE database and sequences associated"

setup(
    name="promgedl",
    version=__version__,

    description=desc,
    url=url,
    author="Antoine Druart",
    author_email="antoine.druart@me.com",

    license="MIT",
    packages=find_packages(),

    entry_point={
        "console_scripts": ["promgedl = promgedl.app:main"],
    }
)