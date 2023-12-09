//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


/*! \file
 *
 * \brief Atmospheric refraction example using an exponential model.
 *
 */


#include "env.hpp"
#include "ray.hpp"

#include "tstModels.hpp"

#include <Engabra>
#include <vector>


namespace
{
	using namespace engabra::g3;

/*
	//! Put current position and tangent values to stream
	std::string
	inline
	nodeStateInfo
		( Vector const & tPrev
		, Vector const & rCurr
		, Vector const & tNext
		, std::size_t const & ndx
		)
	{
		std::ostringstream oss;
		oss
			<< " ndx: " << std::setw(9u) << ndx
			<< " " << "tPrev: " << io::fixed(tPrev, 2u)
			<< " " << "rCurr: " << io::fixed(rCurr, 8u)
			<< " " << "tNext: " << io::fixed(tNext, 2u)
			;
		return oss.str();
	}

	//! Put current node data values to stream
	std::string
	inline
	nodeInfo
		( ray::Node const & node
		, std::size_t const & ndx
		)
	{
		return nodeStateInfo
			(node.thePrevTan, node.theCurrLoc, node.theNextTan, ndx);
	}
*/

} // [anon]


/*! \brief Demonstrate refraction path trace for exponential atmosphere.
 *
 * This is an example used to drive development of various software
 * classes in this project.
 *
 */
int
main
	()
{
	tst::AtmModel const atm(env::sEarth);
	// std::cout << atm.infoString("atm") << std::endl;

	// math/algebra foundation
	using namespace engabra::g3; // for basis vectors, e1,e2,...

	// location on Earth
	double const & groundRad = env::sEarth.theRadGround;
	Vector const stopNear{ groundRad * e3 };

	// initial conditions
	ray::Start const start
		{ ray::Start::from
			( -e3 + .5*e1  // down looking and about 30-deg to the side
			, (groundRad + 9144.)*e3 // about 30k feet altitude
			)
		};

	// ray propgation parms
	constexpr double propStepDist{  10. }; // integration step size
	constexpr double saveStepDist{ 100. }; // save this often

	// path propagation setup
	ray::Propagator const prop{ &atm, propStepDist };
	ray::Path path(start, stopNear, saveStepDist);

	// perform path propagation
	prop.traceNodes(path.theStart, &path);

	// report results
	for (ray::Node const & node : path.theNodes)
	{
		std::cout << node.infoBrief() << std::endl;
	}

	std::cout << path.infoCurvature() << '\n';
}

