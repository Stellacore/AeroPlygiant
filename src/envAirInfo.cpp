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
 * \brief Classes and functions for modeling Air mass properties
 *
 */


#include "envAirInfo.hpp"

#include <Engabra>

#include <fstream>
#include <regex>
#include <sstream>
#include <vector>


namespace aply
{
namespace env
{

	//! Private implementation detail utilities.
	namespace priv
	{
		//! Parse string into a collection of double fields.
		std::vector<double>
		valuesFrom
			( std::string const & str
			)
		{
			std::vector<double> values;
			std::istringstream iss(str);
			// put into temporary until confirming if entire record is okay
			bool allGood{ false };
			std::vector<double> tmpVals;
			tmpVals.reserve(32u); // arbitrary
			double value{};
			while (iss.good())
			{
				value = engabra::g3::null<double>();
				iss >> value;
				if (engabra::g3::isValid(value)) // require all valid values
				{
					tmpVals.emplace_back(value);
					allGood = true;
				}
				else // treat any error as abort
				{
					allGood = false;
					break;
				}
			}
			if (allGood) // if everything okay, then return result.
			{
				values = tmpVals;
			}
			return values;
		}

	} // [priv]


// static
AirInfo
AirInfo :: fromUWyoRecord
	( std::string const & record
		//!< One record line of sounding data
	)
{
	AirInfo info{};

	std::vector<double> const values(priv::valuesFrom(record));
	if (4u < values.size())
	{
		// select values of interest
		double const pres_hPa{ values[0] };
		double const high_J{ values[1] };
		double const temp_C{ values[2] };
	//	double const dewPnt{ values[3] };
		double const relh_pct{ values[4] };

		// convert to standard SI units
			// Assume geopotential height and elevation have
			// essentially equal numeric values (close enough
			// given uncertainty in atomospheric models/data).
		double const high_m{ high_J };
		double const pres_Pa{ 100. * pres_hPa };
		double const temp_K{ 273.15 + temp_C };
		double const relh_f{ .01 * relh_pct };

		// construct return structure
		info = AirInfo{ high_m, temp_K, pres_Pa, relh_f };
	}
	return info;
}

bool
AirInfo :: isValid
	() const
{
	return
		(  engabra::g3::isValid(theHigh)
		&& engabra::g3::isValid(theTemp)
		&& engabra::g3::isValid(thePres)
		&& engabra::g3::isValid(theRelH)
		);
}

double
AirInfo :: temp_C
	() const
{
	return (theTemp - 273.15);
}

double
AirInfo :: pres_mBar
	() const
{
	return (thePres / 100.);
}

double const &
AirInfo :: height
	() const
{
	return theHigh;
}

double
AirInfo :: indexOfRefraction
	() const
{
	return ior::bomford(thePres, theTemp);
}


std::string
AirInfo :: infoBrief
	() const
{
	std::ostringstream oss;
	using engabra::g3::io::fixed;
	oss
		<< fixed(theHigh, 5u, 1u)
		<< " h[m] "
		<< fixed(theTemp, 4u, 1u)
		<< " T[K] "
		<< fixed(thePres, 6u, 0u)
		<< " p[Pa] "
		<< fixed(theRelH, 1u, 3u)
		<< " relH[-] "
		<< fixed(indexOfRefraction(), 1u, 9u)
		<< " IoR[-] "
		;
	return oss.str();
}

	//! \brief Load University WY atmospheric sounding data
	std::map<Height, aply::env::AirInfo>
	airInfoFromUWyoSounding
		( std::ifstream & istrm
		)
	{
		std::map<Height, aply::env::AirInfo> mapHighInfo;
		std::string line;
		while ((! istrm.bad()) && (! istrm.eof()))
		{
			getline(istrm, line);
			// prequalify data lines by skipping those with text description
			static std::regex rxNonDigit("[A-Z]");
			if (std::regex_search(line, rxNonDigit))
			{
				continue;
			}
			// attempt constructing an AirData instance with potential record
			using aply::env::AirInfo;
			AirInfo const info{ AirInfo::fromUWyoRecord(line) };
			if (info.isValid())
			{
				mapHighInfo[info.height()] = info;
			}
		}
		return mapHighInfo;
	}

	//! \brief Load University WY atmospheric sounding data
	std::map<Height, aply::env::AirInfo>
	airInfoFromUWyoSounding
		( std::filesystem::path const & inPath
		)
	{
		std::ifstream ifs(inPath.native());
		return airInfoFromUWyoSounding(ifs);
	}



} // [env]
} // [aply]

