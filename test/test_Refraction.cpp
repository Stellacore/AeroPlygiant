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

#include <Engabra>

/*
#include "libgeo/Refraction.h"

#include "libbase/Tests.h"
#include "libgeo/Atmosphere.h"
#include "libgeo/Ellipsoid.h"
#include "libmath/Ellipsoid.h"
#include "libmath/EllipsoidRay.h"
#include "libmath/math.h"
#include "libmath/Ray3D.h"
#include "libmath/Vector3D.h"
*/

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

/*
	// Null is valid case
	geo::Refraction const null;
	bool const nullIsValid(null.isValid());

	// Make sure there's an infoString
	(void)null.infoString("null");
*/

// ExampleStart
	// quantities typical of remote sensing geometry
	double const angleSensor(engabra::g3::pi / 4.0);
	double const highSensor(2000.);
	double const highGround( 500.);

	// convert to geocentric values for refraction computation
	double const radiusEarth{ aply::env::sEarth.theRadGround };
	double const radiusSensor(radiusEarth + highSensor);
	double const radiusGround(radiusEarth + highGround);

/*
	geo::Refraction refractDown(radiusEarth, radiusSensor, angleSensor);
	double const angleGround(refractDown.angleAt(radiusGround));
// ExampleEnd

	geo::Atmosphere earthAtmosphere(geo::AtmosphereParameters::COESA1976());

	// Compute angle from vertical using a version of Snell's law
	double const angleVerticalGround(std::asin
		( earthAtmosphere.indexOfRefraction(highSensor)
		/ earthAtmosphere.indexOfRefraction(highGround)
		* radiusSensor / radiusEarth
		* std::sin(angleSensor)));

	geo::Refraction refractUp(radiusEarth, radiusGround, angleVerticalGround);
	double const gotVal(refractUp.angleAt(radiusSensor));
	double const expVal(-angleGround);

	// Test displacement

	double const gotDisplacement(refractDown.displacementAt(radiusGround));
	double const expDisplacement
		(radiusGround * angleGround);

	// Test zero refraction for zero inclination angle
	geo::Refraction zeroRefract(radiusEarth, radiusSensor, 0.0);
	double const gotZero(zeroRefract.angleAt(radiusGround));
	double const expZero(0.0);

	// Conditional checking
	if (nullIsValid)
	{
		oss << "failure: nullIsValid" << std::endl;
	}

	if (! math::nearlyEquals(gotVal, expVal))
	{
		oss << "failure of symmetry test:" << std::endl;
		oss << "got: " << gotVal << std::endl;
		oss << "exp: " << expVal << std::endl;
	}

	if (! math::nearlyEquals(gotDisplacement, expDisplacement))
	{
		oss << "failure of displacement test:" << std::endl;
		oss << "got: " << gotDisplacement << std::endl;
		oss << "exp: " << expDisplacement << std::endl;
	}

	if (! math::nearlyEquals(gotZero, expZero))
	{
		oss << "failure of zero inclination angle test:" << std::endl;
		oss << "got: " << gotZero << std::endl;
		oss << "exp: " << expZero << std::endl;
	}
*/

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


