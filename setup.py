#!/usr/bin/env python3

from setuptools import setup, find_packages


setup(
    name="bge_toolkit",
    author="Jackie Goldstein",
    author_email="jigold@broadinstitute.org",
    description="BGE-Toolkit built on top of the Hail Ecosystem.",
    packages=find_packages(),
    package_dir={'bge_toolkit': 'bge_toolkit'},
    python_requires=">=3.9",
    install_requires=['hail', 'pandas>=2,<3', 'numpy', 'typer', 'plotnine', 'sphinx', 'sphinxcontrib-typer', 'sphinx_rtd_theme'],
    entry_points={'console_scripts': ['bge-toolkit = bge_toolkit.cli.__main__:main']},
)
