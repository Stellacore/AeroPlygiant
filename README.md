
# Refraction - repository for optical refraction analysis and simulation


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

