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
 * \brief Demonstrate atmospheric profile from University WY sounding data.
 *
 */


#include "envAirInfo.hpp"
#include "envAirProfile.hpp"

#include <filesystem>
#include <iostream>


namespace
{
	//! \brief Program demoAirSoundingData.cpp main application usage.
	struct Usage
	{
		std::filesystem::path theLoadPath;

		explicit
		Usage
			( int argc
			, char * argv[]
			)
		{
			if (2u == argc)
			{
				theLoadPath = argv[1];
			}
			else
			{
				std::cerr <<
					"\nApplication reads Uni WY atmospheric sounding data"
					"\nand reports index of refraction profile computed"
					"\nfrom those pressure and temperature data"
					"\nThe input data format is that of the UWYO web page"
					"\n(E.g. on web page, select-All, copy, then paste into"
					"\na text file)."
					;
				std::cerr << "\n";
				std::cerr << "Usage: <progName> <UWyoDataPageFile>\n";
				std::cerr << "\n";
			}
		}

		//! True if input path exists
		bool
		isValid
			() const
		{
			return std::filesystem::exists(theLoadPath);
		}

	}; // Usage



} // [anon]


/* \brief Compare computed IoR values between UWyo sounding data and COESA1976.
 *
 * Usage:
 * \arg Program takes one argumen which is the path to a file containing
 * atmospheric sounding data (in "select-all/cut-n-paste" format from the
 * University of Wyoming site:
 * http://weather.uwyo.edu/upperair/sounding.html
 *
 * Code loads these "Sounding" data. It also loads "COESA1976" model data
 * (from aply::env::sAirInfoCoesa1976 in include/envAirInfo.hpp). It then
 * loops over range of heights above ground. At each height it:
 * \arg Intrpolates AirIndex parameters (e.g. Temp/Pres) at height
 * \arg Computes IoR using AirIndex values for each
 * \arg Reports the two IoR values as well as difference between them.
 */
int
main
	( int argc
	, char * argv[]
	)
{
	Usage const use(argc, argv);
	if (! use.isValid())
	{
		return 1;
	}

	//
	// Load raw AirInfo data (e.g. pressure, temperature)
	//

	// Load UWyo atmospheric model
	std::map<aply::env::Height, aply::env::AirInfo> const airMapSounding
		{ aply::env::airInfoFromUWyoSounding(use.theLoadPath) };

	// Load COESA1976 model (from hard coded data)
	std::map<aply::env::Height, aply::env::AirInfo> const airMapCoesa1976
		{ aply::env::sAirInfoCoesa1976 };

	//
	// Wrap data in AirProfile interpolator
	//

	// construct a profile corresponding to sounding data (from arg[1])
	aply::env::AirProfile const profileSounding{ airMapSounding };

	// construct a Standard atmosphere profile for comparison
	aply::env::AirProfile const profileCoesa1976{ airMapCoesa1976 };

	// generate a table of IoR value comparisons
	constexpr double maxHeight{ 15000. };
	constexpr double delHeight{  1000. };
	for (double height{0.} ; height < maxHeight ; height += delHeight)
	{
		double const iorSound{ profileSounding.indexOfRefraction(height) };
		double const iorCoesa{ profileCoesa1976.indexOfRefraction(height) };
		double const iorDiff{ iorSound - iorCoesa };
		using engabra::g3::io::fixed;
		std::cout
			<< "  height: " << fixed(height, 6u, 0u)
			<< "  iorSound: " << fixed(iorSound, 1u, 6u)
			<< "  iorCoesa: " << fixed(iorCoesa, 1u, 6u)
			<< "   iorDiff: " << fixed(iorDiff, 1u, 6u)
			<< '\n';
	}

	return 0;
}

