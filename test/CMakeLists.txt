##
## -- CMake build system description
##

# ===
# === Test Programs
# ===


set(mainProgs

	test_test

	# Algorithm checks
#	test_DiffEq # needs work
	test_Hessian
	test_Sylvester

	# env
	test_IndexVolume

	# ray
	test_nextTangentDir
	test_Path
	test_Propagator
	test_roundTrip

	# geom
	test_Cylinder
	test_Interval

	# Atmospheric refraction code from Stellacore
	test_Atmosphere
	test_Refraction

	)

foreach(mainProg ${mainProgs})

	add_executable(${mainProg} ${mainProg}.cpp)
	target_include_directories(
		${mainProg}
		PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}/../include # public interface
		PRIVATE
			${CMAKE_CURRENT_SOURCE_DIR}  # local test code includes
		)
	target_link_libraries(
		${mainProg}
		PRIVATE
			engabra::engabra
			${aProjLib}  # TODO AeroPlygiantLib::AeroPlygiantLib
		)
	add_test(NAME ${mainProg} COMMAND ${mainProg})

endforeach(mainProg)

