// 
// MIT License
// 
// Copyright (c) 2023 Stellacore Corporation
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// 


/*! \file
 *
 * \brief Atmospheric Refraction Example - round trip ray tracing
 *
 */

#include "tst.hpp"

#include "env.hpp"
#include "ray.hpp"
#include "tstModels.hpp"

#include <Engabra>
#include <vector>


/*! \brief Run ray trace forward and backward - check round trip consistency.
 */
int
main
	()
{
	std::ostringstream oss; // test message string

	tst::AtmModel const atm(env::sEarth);

	// initial conditions
	using namespace engabra::g3;
	Vector const tFwdBeg{ direction(e1 + e3) };
	Vector const rFwdBeg{ env::sEarth.theRadGround * e3 };
	ray::Start const fwdStart{ ray::Start::from(tFwdBeg, rFwdBeg) };

Vector const stopNear{ null<Vector>() };

	// ray tracing parameters
	double const propStepDist{ 0.100 }; // meters
	ray::Propagator const prop{ &atm, propStepDist };

	// trace ray forward
	double const saveStepDist{ 0.100 }; // meters
	ray::Path fwdPath(fwdStart, stopNear, saveStepDist);
	fwdPath.reserve(10u * 1024u);
	prop.tracePath(&fwdPath);

	if (fwdPath.theNodes.empty())
	{
		oss << "Failure of forward path size test\n";
	}
	else
	{
		// for testing, compare forward nodes to those with reverse path (next)

		// get last node in forward pass
		ray::Node const & lastNode = fwdPath.theNodes.back();
		Vector const & tRevBeg = -lastNode.theNextTan;
		Vector const & rRevBeg =  lastNode.theCurrLoc;

		ray::Node const & endNode = fwdPath.theNodes.back();
		Vector const & tExp = tFwdBeg;
		Vector const & rExp = rFwdBeg;
		Vector const & tGot = -endNode.theNextTan; // reverse for test compare
		Vector const & rGot =  endNode.theCurrLoc;;

		// trace ray in reverse direction
		ray::Start const revStart{ ray::Start::from(tRevBeg, rRevBeg) };
		ray::Path revPath(revStart, stopNear, saveStepDist);
		prop.tracePath(&revPath);

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
	}

	if (! oss.str().empty())
	{
		std::cerr << " propStepDist: "
			<< io::fixed(propStepDist, 7u, 6u) << '\n';
		std::cerr << " saveStepDist: "
			<< io::fixed(saveStepDist, 7u, 6u) << '\n';
	}

	return tst::finish(oss);
}

