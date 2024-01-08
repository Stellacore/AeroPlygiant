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

	// quantities typical of remote sensing geometry
	constexpr double angleSensor{ engabra::g3::piQtr };
	constexpr double highSensor{  2000. };
	constexpr double highGround{   500. };

	// convert to geocentric values for refraction computation
	double const radEarth{ aply::env::sEarth.theRadGround };
	double const radSen{ radEarth + highSensor };
	double const radGnd{ radEarth + highGround };

	aply::ray::Refraction const refractDown(radEarth, radSen, angleSensor);
	double const angleGround{ refractDown.thetaAngleAt(radGnd) };

std::cout << "highSensor: " << fixed(highSensor) << '\n';
std::cout << "highGround: " << fixed(highGround) << '\n';
std::cout << "radEarth: " << fixed(radEarth) << '\n';
std::cout << "radSen: " << fixed(radSen) << '\n';
std::cout << "radGnd: " << fixed(radGnd) << '\n';

std::cout << "angleSensor: " << fixed(angleSensor) << '\n';
std::cout << "angleGround: " << fixed(angleGround) << '\n';

// ExampleEnd

	aply::env::Atmosphere earthAtmosphere
		{ aply::env::Atmosphere::COESA1976() };

std::cout << "earthAtmosphere: " << earthAtmosphere.infoContents() << '\n';

	// Compute angle from vertical using a version of Snell's law
	double const atSenIoR{ earthAtmosphere.indexOfRefraction(highSensor) };
	double const atGndIoR{ earthAtmosphere.indexOfRefraction(highGround) };
	double const sinAngAtSen{ std::sin(angleSensor) };
	double const angleVerticalGround
		{ std::asin((atSenIoR/atGndIoR) * (radSen/radEarth) * sinAngAtSen) };

	aply::ray::Refraction const refractUp
		(radEarth, radGnd, angleVerticalGround);
	double const gotVal{ refractUp.thetaAngleAt(radSen) };
	double const expVal{ -angleGround };

	// Test displacement

	double const gotDisplacement{ refractDown.displacementAt(radGnd) };
	double const expDisplacement{ radGnd * angleGround };

	// Test zero refraction for zero inclination angle
	aply::ray::Refraction const zeroRefract(radEarth, radSen, 0.0);
	double const gotZero(zeroRefract.thetaAngleAt(radGnd));
	double const expZero(0.0);

	if (! nearlyEquals(gotVal, expVal))
	{
		oss << "failure of symmetry test:" << std::endl;
		oss << "got: " << fixed(gotVal, 1u, 9u) << std::endl;
		oss << "exp: " << fixed(expVal, 1u, 9u) << std::endl;
	}

	if (! nearlyEquals(gotDisplacement, expDisplacement))
	{
		oss << "failure of displacement test:" << std::endl;
		oss << "got: " << fixed(gotDisplacement, 1u, 9u) << std::endl;
		oss << "exp: " << fixed(expDisplacement, 1u, 9u) << std::endl;
	}

	if (! nearlyEquals(gotZero, expZero))
	{
		oss << "failure of zero inclination angle test:" << std::endl;
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


