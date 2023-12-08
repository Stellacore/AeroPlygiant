//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


/*! \file
 *
 * \brief Atmospheric Refraction Example.
 *
 */


#include "env.hpp"
#include "ray.hpp"

#include <Engabra>
#include <vector>


namespace
{
	using namespace engabra::g3;

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
	env::AtmModel const atm(env::sEarth);

	// initial conditions
	engabra::g3::Vector const tBeg
		{ engabra::g3::direction(engabra::g3::e1 + engabra::g3::e3) };
	engabra::g3::Vector const rBeg
		{ env::sEarth.theRadGround * engabra::g3::e3 };

	// ray tracing parameters
	double const nominalLength{ atm.thickness() };

	// trace ray with various step sizes
	for (double delta{100000.} ; .001 < delta ; delta /=10.)
	{
		ray::Propagator const prop{ &atm, delta };
		std::vector<ray::Node> const fwdNodes
			{ prop.nodePath(tBeg, rBeg, nominalLength) };
		std::cout << " delta: " << engabra::g3::io::fixed(delta, 7u, 6u) << " ";
		std::cout << nodeInfo(fwdNodes.back(), fwdNodes.size());
	}
}

