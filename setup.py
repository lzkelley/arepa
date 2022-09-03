"""
"""

from setuptools import setup

with open('requirements.txt') as inn:
    requirements = inn.read().splitlines()

setup(
    name="arepa",
    author="Luke Zoltan Kelley",
    author_email="lzkelley@northwestern.edu",
    packages=['arepa'],
    install_requires=requirements
)
