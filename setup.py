"""
E2MC2: Python wrapper for ATAT EMC2 cluster-expansion Monte Carlo
"""

from os.path import abspath, dirname
from setuptools import setup, find_packages
import unittest

project_dir = abspath(dirname(__file__))

def unit_tests():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='*_test.py')
    return test_suite

setup(
    name='e2mc2',
    version='0.0.0',
    description='Python wrapper for ATAT EMC2 cluster-expansion Monte Carlo',
    long_description="""
    Python wrapper for ATAT EMC2 cluster-expansion Monte Carlo
    """,
    url="https://github.com/smtg-ucl/e2mc2",
    author="Adam J. Jackson",
    author_email="a.j.jackson@physics.org",
    license='GPL v3',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics'
        ],
    keywords='chemistry dft cluster monte carlo',
    packages=find_packages(),
    package_data={'': ['tests/MgAl2O4.json']},
    install_requires=['pandas'], # Packages go here e.g. numpy, ase
    test_suite='setup.unit_tests',
    entry_points={'console_scripts': [
        'make-traj = e2mc2.tools.make_traj.main',
        'add-atoms = e2mc2.tools.add_atoms.main',
        ]},
    )
