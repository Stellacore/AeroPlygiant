##
## -- CMake build system description
##

cmake_minimum_required(VERSION 3.13)

set(aProjName "AeroPlygiant")
set(aProjVersion "0.1.0") # Set this based on features, updates and patches
set(aProjLib "${aProjName}Lib")

project(
	${aProjName}
	VERSION ${aProjVersion}
	DESCRIPTION "Optical Refraction Analysis and Simulation "
	LANGUAGES CXX
	)

include(CTest)
enable_testing()


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

# ===
# === Documentation
# ===
add_subdirectory(doc) # project reference documentation

# ===
# === Library
# ===

message(aProjLib: ${aProjLib})
# add_library(aProjLib
# 	)

# ===
# === Main Programs
# ===

set(mainProgs

	hello

	demoExpAtmosphere
	demoThickPlate

	)

foreach(mainProg ${mainProgs})

	add_executable (${mainProg} ${mainProg}.cpp)
	target_link_libraries(
		${mainProg}
		PRIVATE
			engabra::engabra
		)

endforeach(mainProg)

# ===
# === Unit tests
# ===

set(testProgs

	test_test

	test_IndexVolume
	test_nextTangentDir
	test_Path
	test_Propagator
	test_roundTrip

	)

foreach(testProg ${testProgs})

	add_executable(${testProg} ${testProg}.cpp)
	target_link_libraries(
		${testProg}
		PRIVATE
			engabra::engabra
		)
	add_test(NAME ${testProg} COMMAND ${testProg})

endforeach(testProg)

