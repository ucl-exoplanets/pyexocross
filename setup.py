#!/usr/bin/env python
from setuptools import find_packages
from distutils.core import setup

packages = find_packages(exclude=('tests', 'doc'))
provides = ['pyexocross', ]


requires = []


install_requires = ['numpy',
                    'scipy',
                    'tabulate',
                    'taurex']

console_scripts = ['pyexocross=pyexocross.run:run_pyexocross']


entry_points = {'console_scripts': console_scripts, }

classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Environment :: Win32 (MS Windows)',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: POSIX :: Linux',
    'Operating System :: Unix',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development :: Libraries',
]

# Handle versioning
version = '0.01'

# with open("README.md", "r") as fh:
#     long_description = fh.read()

setup(name='pyexocross',
    author='Ahmed Faris Al-Refaie',
    author_email='ahmed.al-refaie.12@ucl.ac.uk',
    license="BSD",
    version=version,
    description='Python Exocross',
    classifiers=classifiers,
    packages=packages,
    #long_description=long_description,
    #url='https://github.com/ucl-exoplanets/TauREx3_public/',
    #long_description_content_type="text/markdown",
    #keywords = ['exoplanet','retrieval','taurex','taurex3','atmosphere','atmospheric'],
    #include_package_data=True,
    entry_points=entry_points,
    provides=provides,
    requires=requires,
    install_requires=install_requires,)
    #extras_require={
    #    'Plot':  ["matplotlib"], },
    #)