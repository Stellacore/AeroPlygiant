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

#ifndef aply_env_AirInfo_INCL_
#define aply_env_AirInfo_INCL_

/*! \file
 *
 * \brief Classes and functions for modeling Air mass properties
 *
 */


#include <Engabra>

#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
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
		inline
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

	//! Index of Refraction utilities
	namespace ior
	{
		/*! \brief Optical Index of Refraction value via Bomford 1971.
		 *
		 * Value is computed using Bomford's expression as it is quoted by
		 * Gyer (ref 'gyer1996:AtmRefraction' entry in theory/Papers.bib).
		 */
		inline
		double
		bomford
			( double const & pressurePa
			, double const & temperatureK
			)
		{
			double const mBarPres{ pressurePa / 100. };
			double const refractivity{ .000078831 * (mBarPres / temperatureK) };
			return (1. + refractivity);
		}
	} // [ior]

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

		//! True if all data members are valid
		inline
		bool
		isValid
			() const
		{
			return
				(  engabra::g3::isValid(theHigh)
				&& engabra::g3::isValid(theTemp)
				&& engabra::g3::isValid(thePres)
				&& engabra::g3::isValid(theRelH)
				);
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

		//! Height (elevation) in meters.
		inline
		double const &
		height
			() const
		{
			return theHigh;
		}

		//! \brief Optical Index of Refraction value via Bomford 1971.
		inline
		double
		indexOfRefraction
			() const
		{
			return ior::bomford(thePres, theTemp);
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

	}; // AirInfo

	/*! \brief COESA1976 standard atmosphere model.
	 *
	 * Data (most likely) taken from Gyer's 1996 PE&RS paper on refraction.
	 * Gyer (ref 'gyer1996:AtmRefraction' entry in theory/Papers.bib).
	 */
	std::map<double, AirInfo> const sAirInfoCoesa1976
		{ //           [m]     [K]   [Pa]
		  { -1000., AirInfo{ -1000.0, 294.66, 113930. } }
		, {     0., AirInfo{     0.0, 288.16, 101325. } }
		, {  1000., AirInfo{  1000.0, 281.66,  89876. } }
		, {  2000., AirInfo{  2000.0, 275.16,  79501. } }
		, {  3000., AirInfo{  3000.0, 268.67,  70121. } }
		, {  4000., AirInfo{  4000.0, 262.18,  61660. } }
		, {  5000., AirInfo{  5000.0, 255.69,  54048. } }
		, {  6000., AirInfo{  6000.0, 249.20,  47217. } }
		, {  7000., AirInfo{  7000.0, 242.71,  41105. } }
		, {  8000., AirInfo{  8000.0, 236.23,  35651. } }
		, {  9000., AirInfo{  9000.0, 229.74,  30800. } }
		, { 10000., AirInfo{ 10000.0, 223.26,  26500. } }
		, { 11000., AirInfo{ 11000.0, 216.78,  22700. } }
		, { 12000., AirInfo{ 12000.0, 216.66,  19399. } }
		, { 13000., AirInfo{ 13000.0, 216.66,  16579. } }
		, { 14000., AirInfo{ 14000.0, 216.66,  14170. } }
		, { 15000., AirInfo{ 15000.0, 216.66,  12112. } }
		, { 16000., AirInfo{ 16000.0, 216.66,  10353. } }
		, { 17000., AirInfo{ 17000.0, 216.66,   8850. } }
		, { 18000., AirInfo{ 18000.0, 216.66,   7565. } }
		, { 19000., AirInfo{ 19000.0, 216.66,   6467. } }
		, { 20000., AirInfo{ 20000.0, 216.66,   5529. } }
		, { 21000., AirInfo{ 21000.0, 216.66,   4727. } }
		, { 22000., AirInfo{ 22000.0, 216.66,   4042. } }
		, { 23000., AirInfo{ 23000.0, 216.66,   3456. } }
		, { 24000., AirInfo{ 24000.0, 216.66,   2955. } }
		, { 25000., AirInfo{ 25000.0, 216.66,   2527. } }
		, { 26000., AirInfo{ 26000.0, 219.34,   2163. } }
		};

	//! Alias for std::map key associated with "height above ground"
	using Height = double;

	//! \brief Load University WY atmospheric sounding data
	inline
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
	inline
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


namespace
{
	//! Put aply::env::AirInfo data to standard stream.
	inline
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


#endif // aply_env_AirInfo_INCL_

