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
 * \brief Demonstration: tracing bundle of rays through a thick plate.
 *
 */


#include "tstModels.hpp"

#include "ray.hpp"
//#include "save.hpp"

#include <Engabra>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


//! \brief Functions and data specific to each application (e.g. each demo).
namespace app
{
	using namespace engabra::g3;

	//! Generate a bundle of ray starting elements
	inline
	std::vector<ray::Start>
	rayStarts
		( Vector const & station
		)
	{
		std::vector<ray::Start> starts;
		static std::vector<double> const xVals
			{ -1., -.75, -.50, -.25, .00, .25, .50, .75, 1. };
		static std::vector<double> const yVals
			{ -.4, -.2, .0, .2, .4 };
		static std::vector<double> const zVals
			{ -2. };
		starts.reserve(zVals.size() * yVals.size());
		for (double const & xVal : xVals)
		{
			for (double const & yVal : yVals)
			{
				for (double const & zVal : zVals)
				{
					starts.emplace_back
						(ray::Start::from
							( Vector{ xVal, yVal, zVal }
							, station
							)
						);
				}
			}
		}
		return starts;
	}


	//! \brief Application invokation info
	struct Usage
	{
		std::string theSaveName{};

		//! An instance created from command line argument
		inline
		explicit
		Usage
			( int argc
			, char * argv[]
			)
		{
			if (1 < argc)
			{
				theSaveName = argv[1];
			}
		}

		//! Description of expected invokation
		inline
		std::string
		useMessage
			() const
		{
			std::ostringstream oss;
			oss << '\n';
			oss << "Usage: <progname> <saveFileName>\n";
			return oss.str();
		}

		//! True if all info needed for application is present
		inline
		bool
		isValid
			() const
		{
			return (! theSaveName.empty());
		}

	}; // Usage

} // [app]


/*! \brief Trace bundle of rays through optical flat and save results to file.
 *
 * Writes results to file path provided as command line argument.
 * The plate configuration (normal direction, and thickness, and the
 * various optical indices of refraction (before, inside and after) are
 * hard coded in the Slab constructor call. The configuration of rays
 * starting values is also hardcoded (ref app::rayStarts()).
 *
 */
int
main
	( int argc
	, char * argv[]
	)
{
	app::Usage use(argc, argv);
	if (! use.isValid())
	{
		std::cerr << use.useMessage() << std::endl;
		return 1;
	}

	using namespace engabra::g3;

	tst::Slab const media
		( e3   // 'z' normal direction
		, 4.5  // zBeg
		, 5.5  // zEnd
		, 1.0  // nu below
		, 1.5  // nu inside
		, 1.25  // nu above
		);
	// tst:showMedia(media);

	// configuration
	constexpr double propStepDist{ 1./4096. }; // integration step size
	constexpr double saveStepDist{ 1./128. }; // save this often

	// create tracer
	ray::Propagator const prop{ &media, propStepDist };

	// path specification
	Vector const station { 5., 5., 10. };
	Vector const stopNear{ 5., 5., -5. };

	// starting rays to trace
	std::vector<ray::Start> const starts{ app::rayStarts(station) };

	// trace and report each ray
	std::ofstream ofs(use.theSaveName);
	for (ray::Start const & start : starts)
	{
		// interact with data consumer
		ray::Path path(start, stopNear, saveStepDist);
		prop.tracePath(&path);

		// show path info
		for (ray::Node const & node : path.theNodes)
		{
			ofs << node.infoBrief() << std::endl;
		}
		ofs << "\n\n\n";
		std::cout << "path.size: " << path.theNodes.size() << std::endl;
	}

}

