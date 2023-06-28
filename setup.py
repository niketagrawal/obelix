# -*- coding: utf-8 -*-
#                                                     #
#  __author__ = Adarsh Kalikadien                     #
#  __institution__ = TU Delft                         #
#  __contact__ = a.v.kalikadien@tudelft.nl            #

from setuptools import setup, find_packages

setup(
    name='obelix',
    version='0.1.0',
    packages=find_packages(include=['obelix', 'obelix.*']),
    url='github.com/epics-group/obelix',
    license='GPLv3',
    author='Adarsh Kalikadien',
    author_email='a.v.kalikadien@tudelft.nl',
    description='An automated worklfow for generation & analysis of bidentate ligand containing complexes ',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    package_data={'obelix': ['data']},
    long_description=open('docs/README.md').read(),
)