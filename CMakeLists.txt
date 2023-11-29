##
## -- CMake build system description
##

cmake_minimum_required(VERSION 3.13)

set(aProjName "Refraction")
set(aProjVersion "0.1.0") # Set this based on features, updates and patches

project(
	${aProjName}
	VERSION ${aProjVersion}
	DESCRIPTION "Optical Refraction Analysis and Simulation "
	LANGUAGES CXX
	)


message("### CMAKE_MAJOR_VERSION: " ${CMAKE_MAJOR_VERSION})
message("### CMAKE_MINOR_VERSION: " ${CMAKE_MINOR_VERSION})
message("### CMAKE_PATCH_VERSION: " ${CMAKE_PATCH_VERSION})
message("### CMAKE_TWEAK_VERSION: " ${CMAKE_TWEAK_VERSION})
message("### PROJECT: " ${PROJECT})
message("### PROJECT_VERSION: " ${PROJECT_VERSION})


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(engabra REQUIRED NO_MODULE)
message(engabra Found: ${engabra_FOUND})
message(engabra Version: ${engabra_VERSION})

set(mainProgs

	hello

	)

foreach(mainProg ${mainProgs})

# Executable for build system testing
add_executable (${mainProg} ${mainProg}.cpp)

target_link_libraries(
	${mainProg}
	PRIVATE
		engabra::engabra
	)

endforeach(mainProg)

