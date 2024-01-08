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
\brief This file defines the unit tests for envAtmosphere.
*/

#include "envAtmosphere.hpp"

#include "tst.hpp"

#include <Engabra>

#include <sstream>

namespace
{

//! Check temperature.
void
test0
	( std::ostringstream & oss
	)
{
	// Null is valid case
	aply::env::Atmosphere const null;
	bool const nullIsValid(null.isValid());

	// Make sure there's an infoString
	(void)null.infoString("null");

// ExampleStart
	aply::env::Atmosphere const coesa1976
		{ aply::env::Atmosphere::COESA1976() };

	aply::env::AtmosphereParameters const parms
		{ coesa1976.parametersForHeight(8000.0) };

// ExampleEnd

	using engabra::g3::io::fixed;

	constexpr double checkAtElev{ 0. };
	constexpr double expIndex{ 1.000277 };
	double const gotIndex{ coesa1976.indexOfRefraction(checkAtElev) };
	if (! (0. < gotIndex))
	{
		oss << "Failure of positive indexOfRefraction() test\n";
		oss << "gotIndex: " << fixed(gotIndex) << '\n';
	}
	else
	{
		constexpr double tol{ 0.000001 }; // in the noise for real atm
		if (! engabra::g3::nearlyEquals(gotIndex, expIndex, tol) )
		{
			oss << "Failure of index interpolation test\n";
			oss << "expIndex: " << fixed(expIndex, 1u, 9u) << '\n';
			oss << "gotIndex: " << fixed(gotIndex, 1u, 9u) << '\n';
		}
	}

	double const expTemperature(236.23);

	// Conditional checking
	if (nullIsValid)
	{
		oss << "failure: nullIsValid" << std::endl;
	}

	// Check against tabulated value
	if (! engabra::g3::nearlyEquals(parms.theTemp, expTemperature))
	{
		oss << "failure of temperature test:" << std::endl;
		oss << "got: " << parms.theTemp << std::endl;
		oss << "exp: " << expTemperature << std::endl;
	}
}

/*
//! Check pressure.
std::string geoAtmosphere_test2()
{
	// Use for outputting errors
	std::ostringstream oss;

	geo::AtmosphereParameters const & coesa1976Parameters
		(geo::AtmosphereParameters::COESA1976());
	geo::Atmosphere const coesa1976(coesa1976Parameters);
	double const gotPressure(coesa1976.pressure(15000.0));
	double const expPressure(121.12);

	// Check against tabulated value
	if (! math::nearlyEquals(gotPressure, expPressure))
	{
		oss << "failure of pressure test:" << std::endl;
		oss << "got: " << gotPressure << std::endl;
		oss << "exp: " << expPressure << std::endl;
	}

	return oss.str();
}

//! Check index of refraction.
std::string geoAtmosphere_test3()
{
	// Use for outputting errors
	std::ostringstream oss;

	geo::AtmosphereParameters const & coesa1976Parameters
		(geo::AtmosphereParameters::COESA1976());
	geo::Atmosphere const coesa1976(coesa1976Parameters);

	std::vector<double> gotIndex;

	for (int high(0); high <= 25000; high += 1000)
	{
		gotIndex.push_back(coesa1976.indexOfRefraction(high));
	}

	static double const expIndex[] =
		{ 277.19e-6, 251.55e-6, 227.76e-6, 205.74e-6, 185.40e-6, 166.63e-6
		, 149.36e-6, 133.51e-6, 118.97e-6, 105.68e-6,  93.57e-6,  82.55e-6
		,  70.58e-6,  60.32e-6,  51.56e-6,  44.07e-6,  37.67e-6,  32.20e-6
		,  27.53e-6,  23.53e-6,  20.12e-6,  17.20e-6,  14.71e-6,  12.58e-6
		,  10.75e-6,   9.20e-6};

	double const * itExp = expIndex;
	for (std::vector<double>::const_iterator itGot(gotIndex.begin());
		 itGot != gotIndex.end(); ++itGot, ++itExp)
	{

		// Check against tabulated value
		if (! math::nearlyEquals(*itGot, *itExp + 1.0))
		{
			oss << "failure of index test:" << std::endl;
			oss << "got: " << *itGot << std::endl;
			oss << "exp: " << *itExp + 1.0 << std::endl;
			break;
		}
	}

	return oss.str();
}

//! Check value between given values
std::string geoAtmosphere_test4()
{
	// Use for outputting errors
	std::ostringstream oss;

	geo::AtmosphereParameters const & coesa1976Parameters
		(geo::AtmosphereParameters::COESA1976());
	geo::Atmosphere const coesa1976(coesa1976Parameters);

	double const gotTemperature(coesa1976.temperature(3500.0));
	double const expTemperature((268.67 + 262.18) / 2.0);

	if (! math::nearlyEquals(gotTemperature, expTemperature))
	{
		oss << "failure of interpolation test:" << std::endl;
		oss << "got: " << gotTemperature << std::endl;
		oss << "exp: " << expTemperature << std::endl;
	}

	return oss.str();
}
*/

}

int main(void)
{
	std::ostringstream oss;

	test0(oss);

	//tests.add(geoAtmosphere_test2);
	//tests.add(geoAtmosphere_test3);
	//tests.add(geoAtmosphere_test4);
	return tst::finish(oss);
}
