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

#ifndef aply_env_Air_INCL_
#define aply_env_Air_INCL_

/*! \file
 *
 * \brief Classes and functions for modeling Air mass properties
 *
 */


#include <Engabra>

#include <map>
#include <sstream>
#include <string>
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

	//! \brief Simple container for air property values.
	struct AirInfo
	{
		//! Height [m]
		double theHigh{ engabra::g3::null<double>() };

		//! Temperature [K]
		double theTemp{ engabra::g3::null<double>() };

		//! Pressure [Pa]
		double thePres{ engabra::g3::null<double>() };

		//! Relative Humidity Fraction [-]
		double theRelH{ engabra::g3::null<double>() };


		//! \brief Instance populated with values from UWyo.edu sounding data.
		inline
		static
		AirInfo
		fromUWyoRecord
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

		//! \brief Temperature value in Celsius.
		inline
		double
		temp_C
			() const
		{
			return (theTemp - 273.15);
		}

		//! \brief Pressure value in millBar.
		inline
		double
		pres_mBar
			() const
		{
			return (thePres / 100.);
		}

		//! \brief Short description of values.
		inline
		std::string
		infoBrief
			() const
		{
			std::ostringstream oss;
			using engabra::g3::io::fixed;
			oss
				<< " High[m]: " << fixed(theHigh, 5u, 1u)
				<< " Temp[K]: " << fixed(theTemp, 4u, 1u)
				<< " Pres[Pa]: " << fixed(thePres, 6u, 0u)
				<< " RelH[-]: " << fixed(theRelH, 1u, 3u)
				;
			return oss.str();
		}

	}; // AirInfo


	//! \brief Collection of AirInfo representing an atomspheric profile
	struct AirProfile
	{
		using Height = double;
		std::map<Height, AirInfo> theAirValues;

	}; // Air

} // [env]
} // [aply]


namespace
{
	std::ostream &
	operator<<
		( std::ostream & ostrm
		, aply::env::AirInfo const & airInfo
		)
	{
		ostrm << airInfo.infoBrief();
		return ostrm;
	}

} // [anon]


#endif // aply_env_Air_INCL_

