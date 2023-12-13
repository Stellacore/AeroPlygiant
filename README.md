
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
(in Welsh) or something like "ply-ge-ant" (in English).

### Purpose

AeroPlygiant supports analysis and simulation of basic atmospheric
refraction effects such as those encountered in terrestrial surveying
and airborne (or spaceborne) remote sensing applications.


## Applications

In its current form, AeroPlygiant is primarily a development toolbox
with which specific questions can be investigated by custom coding
something using the available classes. Notwithstanding, some of the
demonstration programs may be generally useful more or less as-is. E.g.

* demoExpAtmosphere - program to simulate refracting ray path from
an airborne sensor platform using nominal (very)simplified exponential
decay model for Earth atmospheric index of refraction.

* demoThickPlate - program with which to evaluate refraction path
through a classic optical "thick plate".


## Features

At the moment features include:

* Gathering of useful notes and list of references in Refraction.lyx document.

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

[Software reference documentation](#Software-Reference-Documentation)
is generated as part of the project build process.


## General Use

TODO - incorporate these steps into example program e.g. demoAeroPlygiant.cpp

In general creating a custom ray path involves the following.

* Define a refractive environment:

	* Derive a class from env::IndexVolume that provides the
	index of refraction values through a volume of space. E.g.

		TODO

	* In general, also derive a class from env::ActiveVolume that
	specifies the region of interest (such that ray propagation
	terminates at boundaries of the active volume). E.g.

		TODO

	* Construct an instance of the refractive media volume. E.g.

		tst::Slab const opticalMedium(....);

* Configure initial ray path boundary conditions:

	* Construct a ray::Start object to define the incident tangent
	direction and an initial point on the ray path. e.g.

		ray::Start const start{ ray::Start::from(tanBeg, locBeg) };

	* Construct a ray::Propagator instance. E.g.

		ray::Propagator const prop{ &opticalMedium, propStepDist };

* Use the ray::Propagator to trace as many rays as desired:

	* Construct a ray::Path instance that will hold the traced
	path information (as a collection of ray::Node instance that
	are saved every delta-length along the path). E.g.

		ray::Path aPath(start, saveDeltaDistance);

	* Use the ray::Propagator to trace as many paths as desired. E.g.

		prop.tracePath(&aPath);

* Retrieve path data:

	* Retrieve desired data from ray::Path archive of ray::Node
	instances. E.g.

		for (ray::Node const & node : aPath.theNodes)
		{
			std::cout << node.infoBrief() << std::endl;
		}


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

Code documentation is generated by Doxygen during the build process and
is available for browsing via e.g.:

	$ <favoriteBrowser> <someBuildDir>/doc/html/index.html

### Installation

Current project configuration does NOT support installation (nor package
building) while CMakeLists.txt still under development. (TBD)

