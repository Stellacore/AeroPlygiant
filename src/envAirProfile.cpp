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
\brief Definitions for aply::env::AirProfile
*/


#include "envAirProfile.hpp"
#include "geomInterval.hpp"

#include <Engabra>

#include <sstream>


namespace aply
{
namespace env
{

AirInfo
AirProfile :: airInfoAtHeight
	( double const & height
	) const
{
	AirInfo info{};
	if (isValid())
	{
		std::map<Height, AirInfo>::const_iterator const nextIter
			{ theAirInfoMap.upper_bound(height) };

		if ( (std::cbegin(theAirInfoMap) != nextIter)
		  && (std::cend(theAirInfoMap) != nextIter)
		   )
		{
			std::map<Height, AirInfo>::const_iterator prevIter{ nextIter };
			--prevIter;

// TODO Move interpolation to a function in AirInfo

			AirInfo const & nextAir = nextIter->second;
			AirInfo const & prevAir = prevIter->second;

			double const & prevH(prevIter->first);
			double const & nextH(nextIter->first);

			using geom::Interval;
			Interval const inverval(prevH, nextH);

			// determine fraction of way into interval
			double const frac{ Interval(prevH, nextH).fracAtValue(height) };

			if ((! (frac < 0.)) && (frac < 1.))
			{
				info.theHigh = height;
				info.theTemp = Interval
					(prevAir.theTemp, nextAir.theTemp).valueAtFrac(frac);
				info.thePres = Interval
					(prevAir.thePres, nextAir.thePres).valueAtFrac(frac);
				info.theRelH = Interval
					(prevAir.theRelH, nextAir.theRelH).valueAtFrac(frac);
			}

		}
	}

	return info;
}
/*
// Atmosphere :: parametersForHeight
{
	AtmosphereParameters parms;

	//! find upper bracketing value
	std::map<double, AtmosphereParameters>::const_iterator
		const nextIter(theParms.upper_bound(height));

	if (nextIter != theParms.end() && nextIter != theParms.begin())
	{
		std::map<double, AtmosphereParameters>::const_iterator
			prevIter(nextIter);
		--prevIter;

		double const & prevHeight(prevIter->first);
		double const & nextHeight(nextIter->first);

		geom::Interval const inverval(prevHeight, nextHeight);

		using geom::Interval;
		double const frac
			{ Interval(prevHeight, nextHeight).fracAtValue(height) };

		parms.theHigh = height;
		parms.theTemp = Interval
			(prevIter->second.theTemp, nextIter->second.theTemp)
			.valueAtFrac(frac);
		parms.thePressure = Interval
			(prevIter->second.thePressure, nextIter->second.thePressure)
			.valueAtFrac(frac);
		parms.theIoR = Interval
			(prevIter->second.theIoR, nextIter->second.theIoR)
			.valueAtFrac(frac);
	}

	return parms;
}
*/



double
AirProfile :: indexOfRefraction
	( double const & height
	) const
{
	double ior{ engabra::g3::null<double>() };
	if (isValid())
	{
		AirInfo info{ airInfoAtHeight(height) };
		ior = info.indexOfRefraction();
	}
	return ior;
}

bool
AirProfile :: isValid
	() const
{
	// need at least two values to be able to interpolate.
	return (1u < theAirInfoMap.size());
}


} // [env]
} // [aply]

