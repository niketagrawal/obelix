# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #

from setuptools import setup, find_packages

setup(
    name='open-bidentate-explorer',
    version='0.1.0',
    packages=find_packages(include=['open-bidentate-explorer', 'open-bidentate-explorer.*']),
    url='github.com/epics-group/open-bidentate-explorer',
    license='GPLv3',
    author='Adarsh Kalikadien',
    author_email='a.v.kalikadien@tudelft.nl',
    description='A Python tool for bias-free automated catalyst structure generation and analysis',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    package_data={'open-bidentate-explorer': ['data/']},
    long_description=open('docs/README.md').read(),
)