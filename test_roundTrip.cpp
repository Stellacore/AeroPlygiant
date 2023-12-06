//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


//! \file Atmospheric Refraction Example - round trip ray tracing


#include "env.hpp"
#include "ray.hpp"

#include <Engabra>
#include <vector>


/*! \brief Run ray trace forward and backward - check round trip consistency.
 */
int
main
	()
{
	int istat{ 1 }; // default to failure

	env::AtmModel const atm(env::sEarth);

	// initial conditions
	using namespace engabra::g3;
	Vector const tFwdBeg{ direction(e1 + e3) };
	Vector const rFwdBeg{ env::sEarth.theRadGround * e3 };

	// ray tracing parameters
	double const delta{ 0.100 }; // meters
	double const nomDist{ atm.thickness() };
	ray::Propagator const prop{ &atm, delta };

	// trace ray forward
	std::vector<ray::Node> fwdNodes;
	fwdNodes.reserve(10u * 1024u); // size determines how much tracing
	prop.traceNodes(tFwdBeg, rFwdBeg, &fwdNodes);

	// get last node in forward pass
	ray::Node const & lastNode = fwdNodes.back();
	Vector const & tRevBeg = -lastNode.theNextTan;
	Vector const & rRevBeg =  lastNode.theCurrLoc;

	// trace ray in reverse direction
	std::vector<ray::Node> revNodes;
	revNodes.reserve(fwdNodes.size());
	prop.traceNodes(tRevBeg, rRevBeg, &revNodes);

	ray::Node const & endNode = revNodes.back();
	Vector const & tExp = tFwdBeg;
	Vector const & rExp = rFwdBeg;
	Vector const & tGot = -endNode.theNextTan; // reverse for test compare
	Vector const & rGot =  endNode.theCurrLoc;;

	std::ostringstream oss;
	double const tol
		{ magnitude(rExp) * std::numeric_limits<double>::epsilon() };
	if (! nearlyEquals(tExp, tGot, tol))
	{
		oss << "Failure of tangent round trip test\n";
		oss << "tExp: "  << tExp << '\n';
		oss << "tGot: "  << tGot << '\n';
	}
	if (! nearlyEquals(rExp, rGot, tol))
	{
		oss << "Failure of location round trip test\n";
		oss << "rExp: "  << rExp << '\n';
		oss << "rGot: "  << rGot << '\n';
	}

	if (oss.str().empty())
	{
		istat = 0;
	}
	else
	{
		std::cerr << " delta: " << io::fixed(delta, 7u, 6u) << '\n';
		std::cerr << oss.str() << '\n';
	}

	return istat;
}

