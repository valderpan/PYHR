#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/03

from setuptools import setup,find_packages


# Use the helper
# h = SetupHelper(initfile="jcvi/__init__.py", readmefile="README.md")
# h.check_version(NAME, majorv=3, minorv=6)
# cmdclass = versioneer.get_cmdclass()
# include_dirs = []
# setup_dir = op.abspath(op.dirname(__file__))
# requirements = [x.strip() for x in open(op.join(setup_dir, "requirements.txt"))]
# h.install_requirements(requires=["cython", "numpy"])


setup(
    name="PYHR",
    version="0.4.4",
    author="Haoran Pan",
    author_email="haoranpan@foxmail.com",
    description="A custom genomics and epigenomics package",
    long_description="A self-use package that integrates genomics and epigenetics",
    long_description_content_type="text",
    url="",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)