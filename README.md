
# AeroPlygiant - Atmospheric refraction simulation and analysis

## Project Info

The name, AeroPlygaint, is a catenation of "aero" and "plygiant". "Aero"
as in pertaining to Earth's atmosphere (as in aerial remote sensing) and
"Plygaint" from Welsh language meaning "refraction". I.e. AeroPlygiant is
a package supporting analysis of atmospheric refraction in the context of
air-borne remote sensing.

## Build Info

First install Engabra software package.

	$ sudo apt-get install ${PACKAGE_PATH}/Engabra-0.1.0-Linux.deb # Install

To build Refraction package

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

