//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


//! \file Atmospheric Refraction Example.


#include "env.hpp"
#include "ray.hpp"

#include <Engabra>
#include <vector>


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
		ray::Propagator const prop{ atm, delta };
		std::vector<ray::Node> const fwdNodes
			{ prop.nodePath(tBeg, rBeg, nominalLength) };
		std::cout << " delta: " << engabra::g3::io::fixed(delta, 7u, 6u) << " ";
		std::cout << ray::nodeInfo(fwdNodes.back(), fwdNodes.size());
	}

/*
return 0;

	// propagate ray forward until maxLength
//	double const delta{ .01 };
	ray::Propagator const prop{ atm, delta };
	std::vector<ray::Node> const fwdNodes
		{ prop.nodePath(tBeg, rBeg, nominalLength) };

	// report results
	std::cout << env::sEarth.infoString("env::sEarth") << '\n';
	for (std::size_t nn{0u} ; nn < fwdNodes.size() ; ++nn)
	{
		std::cout << ray::nodeInfo(fwdNodes[nn], nn);
	}
*/

	/*
	std::cout << "        delta: " << io::fixed(delta) << '\n';
	std::cout << "nominalLength: " << io::fixed(nominalLength) << '\n';
	report(std::cout, fwdNodes.back(), fwdNodes.size());
	*/

}

