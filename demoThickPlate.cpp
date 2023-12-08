//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


/*! \file
 * 
 * \brief Demonstration: tracing bundle of rays through a thick plate.
 *
 */


#include "tstModels.hpp"

#include "ray.hpp"
//#include "save.hpp"

#include <Engabra>

#include <iostream>
#include <sstream>



namespace
{
	using namespace engabra::g3;

} // [anon]


/*! \brief Trace bundle of rays through optical flat and save results to file.
 */
int
main
	()
{
	using namespace engabra::g3;

	tst::Slab const media
		( 4.   // xBeg
		, 6.   // xEnd
		, 1.5  // nu inside
		, 1.0  // nu outside
		);
	// tst:showMedia(media);

	// configuration
	constexpr double propStepDist{ 1./4096. }; // integration step size
	constexpr double saveStepDist{ 1./128. }; // save this often

	// create tracer
	ray::Propagator const prop{ &media, propStepDist };

	// path specification
	Vector const station{ 0., 0., 0. };
	Vector const lookAlong{ 1., 1., 1. };
	ray::Start const start{ ray::Start::from(lookAlong, station) };
	Vector const stopNear{ 10., 10., 10. };

	// interact with data consumer
Vector const & tBeg = start.theTanDir;
Vector const & rBeg = start.thePntLoc;
	ray::Path path(tBeg, rBeg, stopNear, saveStepDist);
	prop.traceNodes(path.theBegTan, path.theBegLoc, &path);

	// show path info
	for (ray::Node const & node : path.theNodes)
	{
		std::cout << node.infoBrief() << std::endl;
	}
	std::cout << "path.size: " << path.theNodes.size() << std::endl;

}

