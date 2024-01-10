##
## -- CMake build system description
##

# ===
# === Demonstration Programs
# ===

set(mainProgs

	hello

	demoIntegrate

	demoAeroPlygiant

	demoAirUWyo
	demoExpAtmosphere
	demoHotRoad
	demoThickPlate

	)


foreach(mainProg ${mainProgs})

	add_executable (${mainProg} ${mainProg}.cpp)
	target_include_directories(
		${mainProg}
		PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}/../include # public interface
		PRIVATE
			${CMAKE_CURRENT_SOURCE_DIR} # module specific includes
			${CMAKE_CURRENT_SOURCE_DIR}/..  # project level code (e.g. example)
		)
	target_link_libraries(
		${mainProg}
		PRIVATE
			engabra::engabra
			${aProjLib}  # TODO AeroPlygiantLib::AeroPlygiantLib
		)

endforeach(mainProg)

