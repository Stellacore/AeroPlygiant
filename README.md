
# AeroPlygiant - Atmospheric refraction simulation and analysis

## Project Info

AeroPlygiant is a C++ development toolset for investigating
general optical refraction behavior associated with arbitrarily complex
three-dimensional (3D) atmospheric refraction conditions.

### Name

The name, AeroPlygaint, is a catenation of "aero" and "plygiant". "Aero"
since pertaining to Earth's atmosphere and "plygaint" from Welsh language
meaning "refraction".

Pronunciation is somewhat up for grabs. Recommended is "air-o" 
followed by ["pluh-g-yant"](https://www.howtopronounce.com/welsh/plygiant)
(in Welsh) or perhaps something like "ply-ge-ant" (in English).

### Purpose

AeroPlygiant supports analysis and simulation of basic atmospheric
refraction effects such as those encountered in airborne (and spaceborne)
remote sensing and terrestrial surveying applications.


## Uses and Applications

In its current form, AeroPlygiant is primarily a development toolbox
with which specific questions can be investigated by custom coding
something easily using the available classes. Notwithstanding, some of the
demonstration programs may be generally useful more or less as-is. E.g.

* demoExpAtmosphere - program to simulate refracting ray path from
an airborne sensor platform using nominal (very)simplified exponential
decay model for Earth atmospheric index of refraction.

* demoThickPlate - program with which to evaluate refraction path
through a classic optical "thick plate".


## Features

At the moment features include:

* Gathering of useful notes and references in Refraction.lyx document.

* Ability (with a bit of coding) to propagate a ray path through a
refractive volume of space in which the index of refraction can vary
arbitrarily in all three diemensions.


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

Reference documentation is then available via

	$ <favoriteBrowser> <someBuildDir>/doc/html/index.html

### Installation

Current project configuration does NOT support installation (nor package
building) while CMakeLists.txt still under development.

