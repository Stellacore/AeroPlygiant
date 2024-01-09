//
// Original Code:
//  Copyright (c) 2007 Stellacore Corporation.
//  Donated to AeroPlygiant Open Source project 2024.
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
\brief Compare computed refraction with Manual of Photogrammetry table data.
*/


#include "tst.hpp"

#include "envPlanet.hpp"
#include "rayRefraction.hpp"

#include <Engabra>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <sstream>
#include <vector>


namespace
{

/* Refraction reference data used in second test function

Source: Manual of Photogrammetry, Forth Edition,
        Editor in Chief Chester C Slama
		American Society of Photogrammetry, Falls Church, Va
		1980, pg. 487

Refraction angle in uRad at 45-deg for ground high 0 km as function of
sensor high - source "3" (3rd Edition of M.o.P.)

H[km]    refraction[uRad]

   .5          4.9
  1.0          9.8
  1.5         14.9
  2.0         19.9
  2.5         25.0
  3.0         30.0
  3.5         35.0
  4.0         39.8
  4.5         44.6
  5.0         49.2
  5.5         53.6
  6.0         57.8
  6.5         61.9
  7.0         65.6
  7.5         69.2
  8.0         72.5
  8.5         75.5
  9.0         78.3
  9.5         80.8
 10.0         83.1

*/

namespace todo
{
} // [anon]


/*! \brief Check integration of Gyer Eqn[12]
 *
 */
std::string
test0
	( std::ostringstream & oss
	)
{
	// Use for outputting errors
	oss.precision(20);

	// for convenience
	using engabra::g3::nearlyEquals;
	using engabra::g3::io::fixed;

// ExampleStart

	// Eample similar to low oblique remote sensing geometry
	constexpr double fwdLookAngle{ engabra::g3::piQtr }; // 45-deg off Nadir
	constexpr double highSensor{  9000. }; // [m] - a bit under 30k'
	constexpr double highGround{     0. }; // [m] - "sea level" for MoP compare
	constexpr double fwdExpRefDevAngle{ .000073300 }; // from MoP table

	// convert to geocentric values for refraction computation
	double const radEarth{ aply::env::sEarth.theRadGround };
	double const radSen{ radEarth + highSensor };
	double const radGnd{ radEarth + highGround };

	aply::ray::Refraction const fwdRefract(fwdLookAngle, radSen, radEarth);
	double const fwdThetaAtEnd{ fwdRefract.thetaAngleAt(radGnd) };

	double const fwdGotRefDevAngle
		{ fwdRefract.angularDeviationFromStart(radGnd, fwdThetaAtEnd) };

	constexpr double tolAngle{ .000005 }; // about 1 arc second

	if (! nearlyEquals(fwdGotRefDevAngle, fwdExpRefDevAngle, tolAngle))
	{
		double const fwdDifRefDevAngle{ fwdGotRefDevAngle - fwdExpRefDevAngle };
		oss << "Failure of forward refraction angle test\n";
		oss << "    highSensor: " << fixed(highSensor, 5u, 3u)
			<< "  [m]\n";
		oss << "    highGround: " << fixed(highGround, 5u, 3u)
			<< "  [m]\n";
		oss << "fwdExpRefDevAngle: " << fixed(fwdExpRefDevAngle, 1u, 6u)
			<< "  From MoP {3rd Ed., pg487}\n";
		oss << "fwdGotRefDevAngle: " << fixed(fwdGotRefDevAngle, 1u, 6u)
			<< "  Using COESA1976 Atmosphere model\n";
		oss << "fwdDifRefDevAngle: " << fixed(fwdDifRefDevAngle, 1u, 6u) << '\n';
	}

// ExampleEnd

constexpr bool showDetail{ true };
if (showDetail)
{
	// Cartesian locations in local polar frame
	using namespace engabra::g3;
	Vector const fwdLocBeg{ radSen * e3 };
	Vector const fwdLocEnd
		{ radGnd * std::sin(fwdThetaAtEnd)
		, 0.
		, radGnd * std::cos(fwdThetaAtEnd)
		};

	// compute angular deviation at the sensor (in local nadir frame)
	static Vector const downDir{ -e3 };
	Vector const fwdLocDel{ fwdLocEnd - fwdLocBeg };
	double const fwdDist{ magnitude(fwdLocDel) };
	BiVector const fwdRefDevAngle3D{ (logG2(downDir * fwdLocDel)).theBiv };
	double const fwdRefDevAngle{ magnitude(fwdRefDevAngle3D) };

	std::cout << '\n';
	std::cout << "=============\n";
	std::cout << "highSensor: " << fixed(highSensor) << '\n';
	std::cout << "highGround: " << fixed(highGround) << '\n';
	std::cout << "radEarth: " << fixed(radEarth) << '\n';
	std::cout << "radSen: " << fixed(radSen) << '\n';
	std::cout << "radGnd: " << fixed(radGnd) << '\n';
	std::cout << "fwdLookAngle: " << fixed(fwdLookAngle, 1u, 9u) << '\n';

	std::cout << '\n';
	std::cout
		<< "fwdThetaAtEnd: "
		<< fixed(fwdThetaAtEnd, 1u, 9u)
		<< '\n';

	std::cout << '\n';
	std::cout << "fwdLocBeg: "
		<< fixed(fwdLocBeg, 7u, 3u)
		<< "  mag: " << fixed(magnitude(fwdLocBeg))
		<< '\n';
	std::cout << "fwdLocEnd: "
		<< fixed(fwdLocEnd, 7u, 3u)
		<< "  mag: " << fixed(magnitude(fwdLocEnd))
		<< '\n';
	std::cout << "fwdLocDel: "
		<< fixed(fwdLocDel, 7u, 3u)
		<< "  mag: " << fixed(magnitude(fwdLocDel))
		<< '\n';

	double const fwdDifRefDevAngle{ fwdGotRefDevAngle - fwdExpRefDevAngle };
std::cout << '\n';
std::cout << "     fwdLookAngle: " << fixed(fwdLookAngle, 1u, 9u) << '\n';
std::cout << "   fwdRefDevAngle: " << fixed(fwdRefDevAngle, 1u, 9u) << '\n';
std::cout << "fwdExpRefDevAngle: " << fixed(fwdExpRefDevAngle, 1u, 9u) << '\n';
std::cout << "fwdGotRefDevAngle: " << fixed(fwdGotRefDevAngle, 1u, 9u) << '\n';
std::cout << "fwdDifRefDevAngle: " << fixed(fwdDifRefDevAngle, 1u, 9u) << '\n';

	std::cout << "=============\n";
	std::cout << '\n';
}

	//
	// Setup ray tracing in reverse
	//

/*
	constexpr double bckLookAngle
		{ engabra::g3::piQtr }; // 45-deg off Nadir

	aply::ray::Refraction const bckRefract(radSen, radEarth, bckLookAngle);
	double const fwdThetaAtEnd{ bckRefract.thetaAngleAt(radGnd) };
*/

	//
	//
	//
	//
	//
	//

	aply::env::Atmosphere earthAtmosphere
		{ aply::env::Atmosphere::COESA1976() };

// std::cout << "earthAtmosphere: " << earthAtmosphere.infoContents() << '\n';

	// Compute angle from vertical using a version of Snell's law
	double const atSenIoR{ earthAtmosphere.indexOfRefraction(highSensor) };
	double const atGndIoR{ earthAtmosphere.indexOfRefraction(highGround) };
	double const sinAngAtSen{ std::sin(fwdLookAngle) };
	double const angleVerticalGround
		{ std::asin((atSenIoR/atGndIoR) * (radSen/radEarth) * sinAngAtSen) };

	aply::ray::Refraction const refractUp
		(angleVerticalGround, radGnd, radEarth);
	double const gotVal{ refractUp.thetaAngleAt(radSen) };
	double const expVal{ -fwdThetaAtEnd };

	// Test displacement

	// Test zero refraction for zero inclination angle
	aply::ray::Refraction const zeroRefract(0., radSen, radEarth);
	double const gotZero(zeroRefract.thetaAngleAt(radGnd));
	double const expZero(0.0);

	if (! nearlyEquals(gotVal, expVal))
	{
		double const difVal{ gotVal - expVal };
		double const relVal{ difVal / expVal };
		oss << "Failure of symmetry test:" << std::endl;
		oss << "got: " << fixed(gotVal, 1u, 9u) << std::endl;
		oss << "exp: " << fixed(expVal, 1u, 9u) << std::endl;
		oss << "dif: " << fixed(difVal, 1u, 18u) << std::endl;
		oss << "rel: " << fixed(relVal, 1u, 18u) << std::endl;
	}

	if (! nearlyEquals(gotZero, expZero))
	{
		oss << "Failure of zero inclination angle test:" << std::endl;
		oss << "got: " << fixed(gotZero, 1u, 9u) << std::endl;
		oss << "exp: " << fixed(expZero, 1u, 9u) << std::endl;
	}

	oss << "\nFailure: restore test0\n";

	return oss.str();
}

}

/*! \brief Unit test for Refraction computation
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);
//	test1(oss); // TODO add from /tmp/dk.tmp/test1.cpp

	return tst::finish(oss);
}


