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

message("adding src library aProjLib: ${aProjLib}")

# set(EXPORT_ALIAS_NAME ${aProjLib}::${aProjLib})
# set(EXPORT_TARGETS_NAME ${aProjLib}Targets)
# set(EXPORT_FILENAME ${EXPORT_TARGETS_NAME}.cmake)
# message("### EXPORT_ALIAS_NAME: " ${EXPORT_ALIAS_NAME})
# message("### EXPORT_TARGETS_NAME: " ${EXPORT_TARGETS_NAME})
# message("### EXPORT_FILENAME: " ${EXPORT_FILENAME})

add_subdirectory(src)

# ===
# === Demonstrations
# ===

add_subdirectory(demo)

# ===
# === Test programs
# ===

add_subdirectory(test)

