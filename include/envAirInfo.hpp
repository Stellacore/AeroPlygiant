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
#include <iostream>
#include <map>
#include <string>
#include <utility>


namespace aply
{
namespace env
{
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
		static
		AirInfo
		fromUWyoRecord
			( std::string const & record
				//!< One record line of sounding data
			);

		//! \brief Interpolate parameters between two samples
		static
		AirInfo
		airInfoInterpolated
			( AirInfo const & beg
			, AirInfo const & end
			, double const & valueAt
			, std::pair<double, double> const & valueBegEnd
			);


		//! True if all data members are valid
		bool
		isValid
			() const;

		//! \brief Temperature value in Celsius.
		double
		temp_C
			() const;

		//! \brief Pressure value in millBar.
		double
		pres_mBar
			() const;

		//! Height (elevation) in meters.
		double const &
		height
			() const;

		//! \brief Optical Index of Refraction value via Bomford 1971.
		double
		indexOfRefraction
			() const;

		//! \brief Short description of values.
		std::string
		infoBrief
			() const;

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
	std::map<Height, aply::env::AirInfo>
	airInfoFromUWyoSounding
		( std::ifstream & istrm
		);

	//! \brief Load University WY atmospheric sounding data
	std::map<Height, aply::env::AirInfo>
	airInfoFromUWyoSounding
		( std::filesystem::path const & inPath
		);



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

