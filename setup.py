"""
tbTAMR --- Genotypic AMR for M. tuberculosis Public Health
"""
from sys import exit, version_info
from setuptools import setup, find_packages
from os import environ
import logging
import tbtamr

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="tbtamr",
    version="1.0.0",
    description="Running TB-Profiler for MDU",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MDU-PHL/tbtamr",
    author="Kristy Horan",
    author_email="kristyhoran15@gmail.com",
    maintainer="Kristy Horan",
    maintainer_email="kristyhoran15@gmail.com",
    python_requires=">=3.10, <4",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    zip_safe=False,
    install_requires=[
                        "pandas",
                        "pytest",
                        "tabulate",
                        "unidecode",
                        "pysam",
                        'requests',
                        "joblib",
                        "pydantic",
                        "tqdm"
                        ],
    test_suite="nose.collector",
    tests_require=["nose", "pytest","psutil"],
    entry_points={
        "console_scripts": [
            "tbtamr=tbtamr.tbtamr:main",
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: Implementation :: CPython",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    package_data={"tbtamr": ["db/*","configs/*"]}
)
