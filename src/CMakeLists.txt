##
## -- CMake build system description
##

# ===
# === Library Components
# ===

set(srcFiles

	envAtmosphereParameters.cpp
	envAtmosphere.cpp

	envAirProfile.cpp
	mathDiffEqSolve.cpp
	rayRefraction.cpp

	)

message("### aProjLib: ${aProjLib}")

add_library(
	${aProjLib}
	STATIC
		${srcFiles}
	)
target_include_directories(
	${aProjLib}
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/../include # public interface
	PRIVATE
		# engabra::engabra # adds literal "engabra::engabra" as local file
		# Perhaps because engabra is not propertly exporting headers??
		/tmpLocal/include/engabra/
	)

