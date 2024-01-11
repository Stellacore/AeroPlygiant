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

#include <Engabra>

#include <iterator>
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
		// get first occurring entry not less than 'height'
		std::map<Height, AirInfo>::const_iterator const nextIter
			{ theAirInfoMap.upper_bound(height) };

		// find the previous entry
		// Note: If nextIter is the first in collection, then the
		// desired "previous" one would be off the end which would
		// represent (extrapolation).
		if ( (std::cbegin(theAirInfoMap) != nextIter) // next is NOT first
		  && (std::cend(theAirInfoMap) != nextIter) // next is in collection
		   )
		{
			// get the entry before this for use in extrapolation
			std::map<Height, AirInfo>::const_iterator
				const prevIter { std::prev(nextIter, 1) };

			// data values to be interpolated
			AirInfo const & prevAir = prevIter->second;
			AirInfo const & nextAir = nextIter->second;

			// end points defining the interpolation
			double const & prevHeight = prevIter->first;
			double const & nextHeight = nextIter->first;

			// interpolate AirInfo values
			info = AirInfo::airInfoInterpolated
				( prevAir, nextAir
				, height
				, std::make_pair(prevHeight, nextHeight)
				);
		}
	}

	return info;
}

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

