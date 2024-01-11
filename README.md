
# AeroPlygiant - Atmospheric refraction simulation and analysis

## Project Info

AeroPlygiant is a C++ development tool set for investigating
general optical refraction behavior associated with arbitrarily complex
three-dimensional (3D) atmospheric refraction conditions.

### Name

The name, AeroPlygaint, is a catenation of "aero" and "plygiant". "Aero"
since pertaining to Earth's atmosphere and "plygaint" from Welsh language
meaning "refraction".

Pronunciation is somewhat up for grabs. Recommended is "air-o" 
followed by ["pluh-g-yant"](https://www.howtopronounce.com/welsh/plygiant)
(in Welsh) or something like "ply-ge-ant" (in English).

### Purpose

AeroPlygiant supports analysis and simulation of basic atmospheric
refraction effects such as those encountered in terrestrial surveying
and airborne (and space) remote sensing applications.


## Applications

In its current form, AeroPlygiant is primarily a development toolbox
with which specific questions can be investigated by custom coding
something using the available classes. Notwithstanding, some of the
demonstration programs may be generally useful more or less as-is. E.g.

### General Case: Finite element propagation

* ./demo/demoAeroPlygiant.cpp - example program that illustrates the
  use of the finite element ray path propagation.

* ./demo/demoHotRoad - program to simulate refraction in a survey sighting
  made along the edge of a road associated with a bubble of low refraction
  index air.

* ./demo/demoExpAtmosphere - program to simulate refracting ray path from
  an airborne sensor platform using nominal (very)simplified exponential
  decay model for Earth atmospheric index of refraction.

* ./demo/demoThickPlate - program with which to evaluate refraction path
  through a classic optical "thick plate".

### Special Case: differential equation solutions

* demo/demoAirSoundingData.cpp - program to compute index of refraction
  using properties of air at different heights within the atmosphere. This
  demonstration computes profiles of IoR values using the both the
  COESA1976 standard atmospheric parameter model and also data downloaded
  from one of the University of WY atmospheric soundings (by reading
  a sounding file specified on command linke). The results of each,
  and the difference between them is reported.

- example program that computes solution
  to atmospheric refraction for ray paths between different height start
  and end points - i.e. classic aerial remote sensing application. The
  solution is based on Gyer's approach (ref theory/Papers.bib).

### Theoretical developments and practice

* demo/demoIntegrate.cpp - example solution to a second order differential
  equation system (as practice and simple analogue for solving the much
  more complex equation system associated with Fermat's theorem).


## Features

At the moment features include:

* Gathering of useful notes and list of references in Refraction.lyx document.

* Ability (with a bit of coding) to propagate a ray path through a
  refractive volume of space in which the index of refraction can vary
  arbitrarily in all three dimensions.

* A simple 4th order Runge Kutta integrator for solving numeric differential
  equation systems.

* Classes and functions for solving differential equation system (Gyer
  model) associated with refraction in the context of classic remote
  sensing applications.


## Resources

Top level project/software description is in README.md (this file).

Technical/math modeling notes are contained in the document
"Refraction.lyx" compatible with the
[Lyx document processor](https://www.lyx.org/).
This project document, along with the associated Papers.bib 
bibliography file provide references to various works on refraction.

[Software reference documentation](#Software-Reference-Documentation)
is generated as part of the project build process.


## General Use

The demonstration program provides a good overview of these steps.

* The demoAeroPlygiant.cpp [example program](demoAeroPlygiant.cpp)
provides a complete example including definition of a custome index
of refraction environment, then tracing and reporting a ray path.

In general producing a custom path trace involves the following steps:

* Define a refractive environment (override aply::env::IndexVolume)
* Specify initial aply::ray::Start geometry (inital condition)
* Use the aply::ray::Propagator to trace path forward
* Access aply::ray::Path data for reporting/output


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
		/repos/AeroPlygiant
	$ cmake --build . --target all -j `nproc`
	$ # -- Run library unit tests
	$ ctest -j `nproc`

### Software Reference Documentation

Code documentation is generated with CMake which uses the
[the Doxygen utility](https://www.doxygen.nl/index.html)
program during the build process to generate HTML documentation. The
resulting documentation can be viewed by e.g.:

	$ <favoriteBrowser> <someBuildDir>/doc/html/index.html

### Installation

Current project configuration does NOT support installation (nor package
building) while CMakeLists.txt still under development. (TBD)

