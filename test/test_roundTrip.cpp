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
#include "envModels.hpp"
#include "ray.hpp"

#include <Engabra>
#include <vector>


/*! \brief Run ray trace forward and backward - check round trip consistency.
 */
int
main
	()
{
	std::ostringstream oss; // test message string

	using namespace aply;
	using namespace engabra::g3;

	env::AtmModel const atm(env::sEarth);
	std::size_t const pathSize{ 8u };

	constexpr double propStepDist{ 1./16. }; // meters
	constexpr double saveStepDist{ 4./16 }; // meters

	// initial conditions
	Vector const tFwdBeg{ direction(e1 + e3) };
	// Pad needs to be larger than the propStepDist.  This is so that
	// the reverse path ray has room to complete it's propgation before
	// being terminated by the (nu=null) value at the bottom of the
	// atmosphere lower boundary.
	double const pad{ 2. * propStepDist };
	Vector const rFwdBeg{ (env::sEarth.theRadGround + pad) * e3 };
	ray::Start const fwdStart{ ray::Start::from(tFwdBeg, rFwdBeg) };

	// ray tracing parameters
	ray::Propagator const prop{ &atm, propStepDist };

	// trace ray forward
	ray::Path fwdPath(fwdStart, saveStepDist);
	fwdPath.reserve(pathSize);
	prop.tracePath(&fwdPath);

/*
std::cout << fwdPath.infoString("fwdPath") << std::endl;
for (ray::Node const & node : fwdPath.theNodes)
{
	std::cout << "    node: " << node.infoBrief() << '\n';
}
std::cout << "fwdPath.size: " << fwdPath.theNodes.size() << std::endl;
*/

	if (fwdPath.theNodes.empty())
	{
		oss << "Failure of forward path size test\n";
	}
	else
	{
		// for testing, compare forward nodes to those with reverse path (next)

		// get last node in forward pass
		ray::Node const & lastNode = fwdPath.theNodes.back();

		ray::Node const revNode{ lastNode.reversed() };

		Vector const & tRevBeg = -lastNode.theNextTan;
		Vector const & rRevBeg =  lastNode.theCurrLoc;

/*
std::cout << '\n';
std::cout << "lastNode: " << lastNode.infoBrief() << '\n';
std::cout << " revNode: " <<  revNode.infoBrief() << '\n';
std::cout << '\n';
*/

		// trace ray in reverse direction
		ray::Start const revStart{ ray::Start::from(tRevBeg, rRevBeg) };
		ray::Path revPath(revStart, saveStepDist);
		revPath.reserve(pathSize);
		prop.tracePath(&revPath);

/*
std::cout << revPath.infoString("revPath") << std::endl;
for (ray::Node const & node : revPath.theNodes)
{
	std::cout << "    node: " << node.infoBrief() << '\n';
}
std::cout << "revPath.size: " << revPath.theNodes.size() << std::endl;
*/

		ray::Node const & revEndNode = revPath.theNodes.back();

		Vector const & tExp = tFwdBeg;
		Vector const & rExp = rFwdBeg;
		Vector const & tGot = -revEndNode.theNextTan;
		Vector const & rGot =  revEndNode.theCurrLoc;;

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

