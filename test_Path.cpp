//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


/*! \file
 * 
 * \brief Simple ray tracing example to test results gathering/save.
 *
 */


#include "tst.hpp"
#include "tstModels.hpp"

#include "env.hpp"
#include "ray.hpp"
#include "save.hpp"

#include <Engabra>

#include <iostream>
#include <sstream>



namespace tst
{
	using namespace engabra::g3;

	//! Evaluate index of refraction on line through media
	void
	showMedia
		( env::IndexVolume const & media
		)
	{
		constexpr double stepDist{ 1./ 4. };
		for (double dist{0.} ; dist < 10. ; dist += stepDist)
		{
			Vector const rCurr{ dist, 5., 5. };
			double const nuCurr{ media.nuValue(rCurr) };
			Vector const gCurr{ media.nuGradient(rCurr, stepDist) };
			std::cout
				<< "  dist: " << io::fixed(dist)
				<< "  nuCurr: " << io::fixed(nuCurr)
				<< "  gCurr: " << io::fixed(gCurr)
				<< std::endl;
		}
	}

} // [tst]


/*! \brief Simple ray trace example with which to check save::Path
 */
int
main
	()
{
	int istat{ 1 }; // default to failure
	std::ostringstream oss;

	using namespace engabra::g3;

	// construct test sphere object
	tst::Sphere const media
		( Vector{ 5., 5.25, 5.125 }
		, 4. // radius
		, 1.1 // nu center
		, 1.0 // nu edge
		);
/*

	tst::Slab const media
		( 4.   // xBeg
		, 6.   // xEnd
		, 1.5  // nu inside
		, 1.0  // nu outside
		);
	// tst:showMedia(media);
*/

	// path specification
	Vector const tBeg{ 5., 5., 5. };
	Vector const rBeg{ 0., 0., 0.  };
	Vector const stopNear{ 10., 10., 10. };

	// configuration
//	constexpr double propStepDist{ 1./4096. }; // integration step size
constexpr double propStepDist{ 1./128. }; // integration step size
	constexpr double saveStepDist{ 1./128. }; // save this often

	// create tracer
	ray::Propagator const prop{ &media, propStepDist };

	// interact with data consumer
	save::Path path(tBeg, rBeg, stopNear, saveStepDist);
	prop.traceNodes(path.theBegTan, path.theBegLoc, &path);

	// show path info
	constexpr bool showIt{ true };
//	constexpr bool showIt{ false };
	if (showIt)
	{
		for (ray::Node const & node : path.theNodes)
		{
			std::cout << node.infoBrief() << std::endl;
		}
		std::cout << "path.size: " << path.theNodes.size() << std::endl;
	}

	// check that path contains data
	if (! (2u < path.size()))  // size of 2 required for test below
	{
		oss << "Failure of path size test\n";
		oss << "path.size: " << path.size() << '\n';
	}
	else
	{
		// check that path exiting the slab is parallel to one entering
		ray::Node const & begNode = path.theNodes.front();
		ray::Node const & endNode = path.theNodes.back();
		Vector const & begDir = begNode.thePrevTan;
		Vector const & endDir = endNode.theNextTan;
		if (! nearlyEquals(begDir, endDir))
		{
			oss << "Failure of begin/end direction test\n";
			oss << "begDir: " << begDir << '\n';
			oss << "endDir: " << endDir << '\n';
		}

		// check that internal ray direction is distinct
		// (assuming block is near center of path)
		std::size_t const midNdx{ path.size() / 2u };
		ray::Node const & midNode = path.theNodes[midNdx];
		Vector const inDir{ .5*(midNode.thePrevTan + midNode.theNextTan) };
		Vector const exDir{ .5*(begNode.thePrevTan + endNode.theNextTan) };
		Spinor const relSpin{ inDir * exDir };
		constexpr double significantDeflection{ .1 };
		double const bivMag{ magnitude(relSpin.theBiv) };
		if (! (significantDeflection < bivMag))
		{
			oss << "Failure of internal significant deflection test\n";
			oss << "inDir: " << inDir << std::endl;
			oss << "exDir: " << exDir << std::endl;
			oss << "relSpin: " << relSpin << std::endl;
			oss << " bivTol: " << significantDeflection << std::endl;
			oss << " bivMag: " << bivMag << std::endl;
		}
	}

	return tst::finish(oss.str());
}

