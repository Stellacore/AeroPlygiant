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
 * \brief Unit test for class geom::Cylinder
 *
 */


#include "tst.hpp"

#include "geomCylinder.hpp"

#include <sstream>


namespace
{
	//! Check .... TODO
	void
	test0
		( std::ostringstream & oss
		)
	{
		// [DoxyExample00]
		using namespace engabra::g3;
		Vector const begAxis{ zero<Vector>() };
		Vector const dirAxis{ e3 };
		double const rad{ 5. };
		double const len{ 17. };

		// construct cylinder with above parameters
		using namespace aply;
		geom::Cylinder const cylinder(begAxis, dirAxis, len, rad);

		// decompose location into values relative to cylinder
		Vector const aLoc{ 2., 0.00, 11. };
		double const gotRadialDist{ cylinder.distanceFromAxis(aLoc) };
		double const gotRadialFrac{ cylinder.fractionFromAxis(aLoc) };
		double const gotLengthDist{ cylinder.distanceAlongAxis(aLoc) };
		double const gotLengthFrac{ cylinder.fractionAlongAxis(aLoc) };
		double const expRadialDist{ std::sqrt(2.*2. + 0.00*0.00) };
		double const expRadialFrac{ expRadialDist/5. };
		double const expLengthDist{ 11. };
		double const expLengthFrac{ 11./17. };
		// [DoxyExample00]

		using tst::checkGotExp;
		checkGotExp(oss, gotRadialDist, expRadialDist, "RadialDist");
		checkGotExp(oss, gotRadialFrac, expRadialFrac, "RadialFrac");
		checkGotExp(oss, gotLengthDist, expLengthDist, "LengthDist");
		checkGotExp(oss, gotLengthFrac, expLengthFrac, "LengthFrac");
	}
}


/*! \brief Unit test for geom::Cylinder
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);

	return tst::finish(oss);
}

