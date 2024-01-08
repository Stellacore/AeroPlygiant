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
\brief Definitions for aply::env::Atmosphere
*/


#include "envAtmosphere.hpp"
#include "geomInterval.hpp"

#include <sstream>


namespace
{

	using namespace aply::env;
	AtmosphereParameters const coesa1976[] = {
		  //                        [m]     [K]   [mBar]    [unitless]
		  AtmosphereParameters{ -1000.0, 294.66, 1139.30, 1.+304.80e-6 }
		, AtmosphereParameters{     0.0, 288.16, 1013.25, 1.+277.19e-6 }
		, AtmosphereParameters{  1000.0, 281.66,  898.76, 1.+251.55e-6 }
		, AtmosphereParameters{  2000.0, 275.16,  795.01, 1.+227.76e-6 }
		, AtmosphereParameters{  3000.0, 268.67,  701.21, 1.+205.74e-6 }
		, AtmosphereParameters{  4000.0, 262.18,  616.60, 1.+185.40e-6 }
		, AtmosphereParameters{  5000.0, 255.69,  540.48, 1.+166.63e-6 }
		, AtmosphereParameters{  6000.0, 249.20,  472.17, 1.+149.36e-6 }
		, AtmosphereParameters{  7000.0, 242.71,  411.05, 1.+133.51e-6 }
		, AtmosphereParameters{  8000.0, 236.23,  356.51, 1.+118.97e-6 }
		, AtmosphereParameters{  9000.0, 229.74,  308.00, 1.+105.68e-6 }
		, AtmosphereParameters{ 10000.0, 223.26,  265.00, 1.+ 93.57e-6 }
		, AtmosphereParameters{ 11000.0, 216.78,  227.00, 1.+ 82.55e-6 }
		, AtmosphereParameters{ 12000.0, 216.66,  193.99, 1.+ 70.58e-6 }
		, AtmosphereParameters{ 13000.0, 216.66,  165.79, 1.+ 60.32e-6 }
		, AtmosphereParameters{ 14000.0, 216.66,  141.70, 1.+ 51.56e-6 }
		, AtmosphereParameters{ 15000.0, 216.66,  121.12, 1.+ 44.07e-6 }
		, AtmosphereParameters{ 16000.0, 216.66,  103.53, 1.+ 37.67e-6 }
		, AtmosphereParameters{ 17000.0, 216.66,   88.50, 1.+ 32.20e-6 }
		, AtmosphereParameters{ 18000.0, 216.66,   75.65, 1.+ 27.53e-6 }
		, AtmosphereParameters{ 19000.0, 216.66,   64.67, 1.+ 23.53e-6 }
		, AtmosphereParameters{ 20000.0, 216.66,   55.29, 1.+ 20.12e-6 }
		, AtmosphereParameters{ 21000.0, 216.66,   47.27, 1.+ 17.20e-6 }
		, AtmosphereParameters{ 22000.0, 216.66,   40.42, 1.+ 14.71e-6 }
		, AtmosphereParameters{ 23000.0, 216.66,   34.56, 1.+ 12.58e-6 }
		, AtmosphereParameters{ 24000.0, 216.66,   29.55, 1.+ 10.75e-6 }
		, AtmosphereParameters{ 25000.0, 216.66,   25.27, 1.+  9.20e-6 }
		, AtmosphereParameters{ 26000.0, 219.34,   21.63, 1.+  7.77e-6 }
	};

	size_t const coesaSize(sizeof(coesa1976) / sizeof(coesa1976[0]));

} // [anon]

namespace aply
{
namespace env
{


Atmosphere
Atmosphere :: COESA1976
	()
{
	Atmosphere res;

	/*
	for (AtmosphereParameters const * iter(coesa1976)
		; iter != (coesa1976 + coesaSize); +++iter)
	*/
	for (size_t index(0); index < coesaSize; ++index)
	{
		AtmosphereParameters const & parms(coesa1976[index]);
		res.theParms[parms.theHigh] = parms;
	}

	return res;
}


Atmosphere :: Atmosphere
	()
	: theParms()
{
}


bool
Atmosphere :: isValid() const
{
	return 1 < theParms.size();
}


AtmosphereParameters
Atmosphere :: parametersForHeight
	( double const & height
	) const
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


double
Atmosphere :: indexOfRefraction
	( double const & height
	) const
{
	AtmosphereParameters parmsAtHeight{ parametersForHeight(height) };
	return parmsAtHeight.theIoR;
}


std::string
Atmosphere :: infoString
	( std::string const & title
	, std::string const & // fmt
	) const
{
	std::ostringstream oss;

	if (!title.empty())
	{
		oss << title << std::endl;
	}

	oss << "Size: " << theParms.size();

	return oss.str();
}

} // [env]
} // [aply]

