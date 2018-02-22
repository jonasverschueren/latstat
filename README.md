# latstat [![Build Status](https://travis-ci.org/jonasverschueren/latstat.svg?branch=master)](https://travis-ci.org/jonasverschueren/latstat)

Calculate phonon dispersion curves from molecular dynamics potential files.

## Requirements
A C++ compiler with C++11 support and a version of CMake no less than 3.4.2. Earlier versions of CMake 3.x.x are expected to work as well but have not been tested. 

## Download & get started
- Clone this git repository and additional software (googletest and cubature) should download during the build process. 
- Build the package in a new build folder recommended and execute cmake ..
- compile : make
- run the included test suite: ./test/test_suite.exe

## Usage
The software comes with an example main ex_main.cpp which shows how to call the EvaluatePhononFrequencySquared function. 
Available crystal structures: Square mesh, simple cubic, body-centred cubic, face-centred cubic. 
Available potentials: Harmonic and Cosinusoidal 1D interactions, Finnis-Sinclair potentials (see pot/W87.eam.fs), any mono-atomic EamSetfl potential file as used in LAMMPS. 
The code should be easy to extend and adapt to taste.

## Third party software
This software package is distributed with a copy of the Eigen Matrix Library (eigen.tuxfamily.org) to avoid version clashing. In the build process the googletest.git and cubature.git repositories are cloned or kept up-to-date. The latter is included to facilitate integrations over Brillouin zones, a common occurence in phonon calculations.
