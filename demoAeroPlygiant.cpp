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
 * \brief Demonstrate a basic AeroPlygiant ray path solution.
 *
 */



#include "envActiveVolume.hpp"
#include "envIndexVolume.hpp"
#include "rayPath.hpp"
#include "rayPropagator.hpp"
#include "rayStart.hpp"

#include <iostream>


namespace
{
	//! \brief Example of optical refractive medium 3D region
	struct Media : public aply::env::IndexVolume
	{
		explicit
		Media
			( std::shared_ptr<aply::env::ActiveVolume> const & ptVolume
			)
			: IndexVolume(ptVolume)
		{ }

		//! Create a double convex lens in "air" (with decentered surfaces)
		inline
		double
		nuValue
			( engabra::g3::Vector const & rVec
			) const
		{
			using namespace engabra::g3;
			// set index of refraction value based on optical system geometry
			double nu{ 1. }; // air like
			constexpr double r1Sq{ sq(11.) };
			static Vector const c1{ -10., 0., 0. };
			constexpr double r2Sq{ sq(21.) };
			static Vector const c2{  20., 0., 1. };
			if ((magSq(rVec - c1) < r1Sq) && (magSq(rVec - c2) < r2Sq))
			{
				// inside the "lens" set IoR to be glass-like
				nu = 1.5;
			}
			return nu;
		}

	}; // Media
}


/*! \brief Unit test for CN
 */
int
main
	()
{
	// define working volume
	std::shared_ptr<aply::env::ActiveVolume> const ptVolume
		{ std::make_shared<aply::env::ActiveBox>
			( engabra::g3::Vector{ -5., -10., -10. }
			, engabra::g3::Vector{  5.,  10.,  10. }
			)
		};
	// use volume to "clip" media
	Media const media(ptVolume);

	// specify initial conditions
	using namespace engabra::g3;
	Vector const tBeg{ direction(Vector{ 1., .2, .3 }) };
	Vector const rBeg{ -5., 0., 0. };
	Vector const approxEndLoc{ 10., 0., 0. };
	aply::ray::Start const start{ aply::ray::Start::from(tBeg, rBeg) };

	// propagate ray forward
	constexpr double propStepDist{ 1./1024. };
	constexpr double saveStepDist{ 1./16. };

	// configure propagor and trace path
	aply::ray::Propagator const prop{ &media, propStepDist };
	aply::ray::Path path(start, saveStepDist, approxEndLoc);
	prop.tracePath(&path);

	// display ray path (saved nodes)
	for (aply::ray::Node const & node : path.theNodes)
	{
		std::cout << node.infoBrief() << std::endl;
	}
	std::cout << "# path.size: " << path.theNodes.size() << std::endl;

	return 0;
}

