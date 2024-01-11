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
#include "envAtmosphere.hpp"

#include <filesystem>
#include <iostream>


namespace
{
	//! demoAir application usage
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
	// Wrap data in Atmosphere properties interpolation class
	//

	aply::env::AirProfile const profileSounding{ airMapSounding };
	aply::env::AirProfile const profileCoesa1976{ airMapCoesa1976 };

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

return 1;


// TODO do something with this
	// Create COESA standard model
	aply::env::Atmosphere const coesa1976
		{ aply::env::Atmosphere::COESA1976() };

	std::cout << "# loaded from: " << use.theLoadPath << '\n';
	for (std::map<aply::env::Height, aply::env::AirInfo>::value_type
		const & pairHighInfo : airMapSounding)
	{
		aply::env::AirInfo const & info = pairHighInfo.second;
		double const uwyoIoR{ info.indexOfRefraction() };
		double const height{ info.height() };
		double const coesaIoR{ coesa1976.indexOfRefraction(height) };
		double const diff{ uwyoIoR - coesaIoR };
		using engabra::g3::io::fixed;
//		std::cout << info << '\n';
		std::cout
			<< "  height: " << fixed(height, 6u, 0u)
			<< "  uWyoIor: " << fixed(uwyoIoR, 1u, 6u)
			<< "  coesa: " << fixed(coesaIoR, 1u, 6u)
			<< "  diff: " << fixed(diff, 1u, 6u)
			<< '\n';
	}
}

