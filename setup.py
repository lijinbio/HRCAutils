#!/usr/bin/env python3

import os
import io
import re
import ast

from setuptools import setup, find_packages

DEPENDENCIES=['click>=8.1', 'importlib-metadata', 'pathlib']
EXCLUDE_FROM_PACKAGES=["contrib", "docs", "tests*"]
CURDIR=os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(CURDIR, "README.md"), "r", encoding="utf-8") as readme_file:
    readme=readme_file.read()

def get_version():
    main_file=os.path.join(CURDIR, "src", "HRCAutils", "__init__.py")
    _version_re=re.compile(r"__version__\s+=\s+(?P<version>.*)")
    with open(main_file, "r", encoding="utf8") as f:
        match=_version_re.search(f.read())
        version=match.group("version") if match is not None else '"unknown"'
    return str(ast.literal_eval(version))

setup(
    name="HRCAutils",
    author="Jin Li",
    author_email="lijin.abc@gmail.com",
    python_requires=">=3.10",
    description="Utilities of HRCA v1.0.",
    install_requires=DEPENDENCIES,
    packages=find_packages(
        where='src',
        # include=['pkg*'],
        exclude=EXCLUDE_FROM_PACKAGES
        ),
    package_dir={'': 'src'},
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords=[],
    scripts=[],
    setup_requires=[],
    entry_points={
        "console_scripts": [
            "HRCAutils=HRCAutils.HRCAutils:main", # main package
            "exeBash=exeBash.exeBash:exeBash",
            "exePython=exePython.exePython:exePython",
            "scarchesh5ad2refbyscanvi=scarchesh5ad2refbyscanvi.scarchesh5ad2refbyscanvi:main",
            "scarches2queryannobyscanvi=scarches2queryannobyscanvi.scarches2queryannobyscanvi:main",
            ],
        },
    url="https://github.com/lijinbio/HRCAutils",
    version=get_version(),
    zip_safe=False,
    license="MIT license",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)
