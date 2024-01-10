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
#include "envAtmosphere.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>


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


	using Height = double;

	//! \brief Load University WY atmospheric sounding data
	std::map<Height, aply::env::AirInfo>
	airMapUWyoSoundingFrom
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
	airMapUWyoSoundingFrom
		( std::filesystem::path const & inPath
		)
	{
		std::ifstream ifs(inPath.native());
		return airMapUWyoSoundingFrom(ifs);
	}


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

	// Load UWyo atmospheric model
	std::map<Height, aply::env::AirInfo> const airMapSounding
		{ airMapUWyoSoundingFrom(use.theLoadPath) };

	// Load COESA1976 model (from hard coded data)
	std::map<Height, aply::env::AirInfo> const airMapCoesa1976
		{ aply::env::sAirInfoCoesa1976 };

// TODO do something with this
	// Create COESA standard model
	aply::env::Atmosphere const coesa1976
		{ aply::env::Atmosphere::COESA1976() };

	std::cout << "# loaded from: " << use.theLoadPath << '\n';
	for (std::map<Height, aply::env::AirInfo>::value_type
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

