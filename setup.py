# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# Software de Observaciones Sintéticas" S.O.S.
# Setup tools
#
# SPA group at INAOE, @ 24 August 2020
# Latest Revision: 24 Aug 2020, 17:39 GMT
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sospy-pkg-INAOE", # Replace with your own username
    version="0.0.1",
    author="SPA group",
    author_email="mbecerrilt92@gmail.com",
    description="Software de Observaciones Sintéticas sospy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",    # GitHub link
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
