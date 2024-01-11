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
 * \brief Simple ray tracing example to test results gathering/save.
 *
 */


#include "tst.hpp"

#include "env.hpp"
#include "ray.hpp"

#include "example/indexModel.hpp"

#include <Engabra>

#include <iostream>
#include <sstream>



namespace tst
{
	using namespace engabra::g3;

	//! Evaluate index of refraction on line through media
	void
	showMedia
		( aply::env::IndexVolume const & media
		)
	{
		constexpr double stepDist{ 1./ 4. };
		for (double dist{0.} ; dist < 10. ; dist += stepDist)
		{
			Vector const rCurr{ dist, 5., 5. };
			double const nuCurr{ media.nuValue(rCurr) };
			Vector const gCurr{ media.nuGradient(rCurr, stepDist) };
			std::cout
				<< "  dist: " << io::fixed(dist)
				<< "  nuCurr: " << io::fixed(nuCurr)
				<< "  gCurr: " << io::fixed(gCurr)
				<< std::endl;
		}
	}

} // [tst]


namespace
{

	//! Check propagation through a classic thick plate
	std::string
	test0
		()
	{
		std::ostringstream oss;

		using namespace aply;
		using namespace engabra::g3;

		std::shared_ptr<env::ActiveVolume> const ptVolume
			{ std::make_shared<env::ActiveBox>
				(zero<Vector>(), Vector{10., 10., 10.})
			};
		aply::env::index::Slab const media
			( e1   // 'x' direction
			, 4.   // xBeg
			, 6.   // xEnd
			, 1.0  // nu before
			, 1.5  // nu inside
			, 1.0  // nu after
			, ptVolume // a bounding box
			);
		// tst:showMedia(media);

		// path specification
		Vector const tBeg{ 5., 5., 5. };
		Vector const rBeg{ 0., 0., 0.  };
		Vector const approxEndLoc{ 10., 10., 10. };
		ray::Start const start{ ray::Start::from(tBeg, rBeg) };

		// configuration
		constexpr double propStepDist{ 1./128. }; // integration step size
		constexpr double saveStepDist{ 1./128. }; // save this often

		// create tracer
		ray::Propagator const prop{ &media, propStepDist };

		// interact with data consumer
		ray::Path path(start, saveStepDist, approxEndLoc);
		prop.tracePath(&path);

		/*
		std::cout
			<< "\nCompleted path:\n"
			<< ray::PathView{&path}.infoCurvature()
			<< '\n';
		*/

		// show path info
	//	constexpr bool showIt{ true };
		constexpr bool showIt{ false };
		if (showIt)
		{
			for (ray::Node const & node : path.theNodes)
			{
				std::cout << node.infoBrief() << std::endl;
			}
			std::cout << "path.size: " << path.theNodes.size() << std::endl;
		}

		// check that path contains data
		if (! (2u < path.size()))  // size of 2 required for test below
		{
			oss << "Failure of path size test\n";
			oss << "path.size: " << path.size() << '\n';
		}
		else
		{
			// check that path exiting the slab is parallel to one entering
			ray::Node const & begNode = path.theNodes.front();
			ray::Node const & endNode = path.theNodes.back();
			Vector const & begDir = begNode.thePrevTan;
			Vector const & endDir = endNode.theNextTan;
			if (! nearlyEquals(begDir, endDir))
			{
				oss << "Failure of begin/end direction test\n";
				oss << "begDir: " << begDir << '\n';
				oss << "endDir: " << endDir << '\n';
			}

			// check that internal ray direction is distinct
			// (assuming block is near center of path)
			std::size_t const midNdx{ path.size() / 2u };
			ray::Node const & midNode = path.theNodes[midNdx];
			Vector const inDir{ .5*(midNode.thePrevTan + midNode.theNextTan) };
			Vector const exDir{ .5*(begNode.thePrevTan + endNode.theNextTan) };
			Spinor const relSpin{ inDir * exDir };
			constexpr double significantDeflection{ .1 };
			double const bivMag{ magnitude(relSpin.theBiv) };
			if (! (significantDeflection < bivMag))
			{
				oss << "Failure of internal significant deflection test\n";
				oss << "inDir: " << inDir << std::endl;
				oss << "exDir: " << exDir << std::endl;
				oss << "relSpin: " << relSpin << std::endl;
				oss << " bivTol: " << significantDeflection << std::endl;
				oss << " bivMag: " << bivMag << std::endl;
			}
		}

		return oss.str();
	}

} // [anon]

/*! \brief Simple ray trace example with which to check ray::Path
 */
int
main
	()
{
	std::ostringstream oss;

	oss << test0();

	return tst::finish(oss.str());
}

