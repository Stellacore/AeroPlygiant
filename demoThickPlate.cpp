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

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>



namespace
{
	using namespace engabra::g3;

	//! Generate a bundle of ray starting elements
	std::vector<ray::Start>
	rayStarts
		( Vector const & station
		)
	{
		std::vector<ray::Start> starts;
		static double const xVal{ .5 };
		static std::vector<double> const yVals{ .00, .25, .50, .75, 1. };
		static std::vector<double> const zVals{ .00, .25, .50, .75, 1. };
		starts.reserve(zVals.size() * yVals.size());
		for (double const & yVal : yVals)
		{
			for (double const & zVal : zVals)
			{
				starts.emplace_back
					(ray::Start::from(Vector{ xVal, yVal, zVal }, station));
			}
		}
		return starts;
	}

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
	Vector const station{  0.,  0.,  0. };
	Vector const stopNear{ 10., 10., 10. };

	// starting rays to trace
	std::vector<ray::Start> const starts{ rayStarts(station) };

	// trace and report each ray
std::ofstream ofs("/tmp/dk.tmp/foo");
	for (ray::Start const & start : starts)
	{
		// interact with data consumer
		ray::Path path(start, stopNear, saveStepDist);
		prop.traceNodes(path.theStart, &path);

		// show path info
		for (ray::Node const & node : path.theNodes)
		{
			ofs << node.infoBrief() << std::endl;
		}
		ofs << "\n\n\n";
		std::cout << "path.size: " << path.theNodes.size() << std::endl;
	}

}

