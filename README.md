
# AeroPlygiant - Atmospheric refraction simulation and analysis


## Project Info

The name, AeroPlygaint, is a catenation of "aero" and "plygiant". "Aero"
(air-o) in pertaining to Earth's atmosphere and "plygaint" (plug-yant)
from Welsh language meaning "refraction".

AeroPlygiant supports analysis and simulation of basic atmospheric
refraction effects that are encountered in airborne (and spaceborne)
remote sensing and terrestrial surveying applications.


## Uses and Applications

In its current form, AeroPlygiant is primarily a development toolbox
with which specific questions can be investigated by custom coding
something to answer the question at hand. Notwithstanding, some of the
demonstration programs may be generally useful as command line utility
applications. E.g.

* (TODO)


## Resources

Top level project/software description is in README.md (this file).

Technical/math modeling notes are contained in the document
"Refraction.lyx" compatible with the
[Lyx document processor](https://www.lyx.org/).
This project document, along with the associated \ref Papers.bib 
bibliography file provide references to various works on refraction.


## Project Development Environment

The software is expressed in C++ (standard 17) and should be fairly
portable. It is known to build with:

* [g++](https://gcc.gnu.org/) (GCC) 12.2.0 - Gnu Compiler Collection
for code compilation, assembly and linking.

* [CMake](https://cmake.org/) - Build system generator. Used to
create systems for project build (CMake), testing(CTest), installation,
and packaging(CPack)

Build dependencies:

* [Engabra - the ENGineering Geometric AlegBRA pacakge](
	https://github.com/Stellacore/Engabra/)
	- Geometric algebra library used extensively for math/geometry
	computations.

Development dependencies

* [Doxygen](https://www.doxygen.nl/) - Reference document generation.

* [Lyx](https://www.lyx.org/) - Document processor application
that utilizes Tex/LaTex to produce print-ready quality documents in
various formats.

* [Ipe](https://ipe.otfried.org/) - Extensible drawing editor (used
for producing figures for Lyx documents.


## Build Info

First install the required Engabra software package.

	$ sudo apt-get install ${PACKAGE_PATH}/Engabra-0.1.0-Linux.deb # Install

To build Refraction package. For example (replace specific directory
names and included options as appropriate for local system).

	$ # -- Compile everything (including doxygen documentation)
	$ mkdir <someBuildDir> && cd <someBuildDir>
	$ cmake  \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_INSTALL_PREFIX=/tmpLocal/ \
		-DCMAKE_PREFIX_PATH=/tmpLocal/ \
		/repos/Refraction
	$ cmake --build . --target all -j `nproc`
	$ # -- Run library unit tests
	$ ctest -j `nproc`

### Installation

Current project configuration does NOT support installation (nor package
building) while CMakeLists.txt still under development.

