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
 * \brief Demonstration example of basic AeroPlygiant use.
 *
 */


#include "envActiveVolume.hpp" // Defines limits of ray propagation
#include "envIndexVolume.hpp" // Base class for expressing IoR media
#include "rayPath.hpp" // Data structure for storing propagated ray path
#include "rayPropagator.hpp" // Functions to compute ray path
#include "rayStart.hpp" // Initial conditions (tangent dir and location)

#include <functional>
#include <iostream>


namespace
{
	// [DoxyExampleMedia]
	/*! \brief Example of optical refractive medium 3D region.
	 *
	 * Overload the (abstract) baseclass to implement any arbitrary
	 * Index of Refraction (IoR) scalar field.
	 *
	 * Note that the baseclass accepts a shared pointer to an
	 * env::ActiveVolume instance, which is used to "clip" the media
	 * to this volume (ray propagation computations will stop when
	 * encountering the edge of the active volume).
	 */
	struct Media : public aply::env::IndexVolume
	{
		//! Attach env::ActiveVolume to define boundaries of the IoR field.
		explicit
		Media
			( std::shared_ptr<aply::env::ActiveVolume> const & ptVolume
			)
			: IndexVolume(ptVolume)
		{ }

		/*! \brief Specify the IoR scalar field specific to problem at hand.
		 *
		 * For this example, model a double convex lens in air. The "lens"
		 * is formed computationally by considering the intersection of
		 * two spheres. When inside this intersection, an IoR value of
		 * 1.500 is returned (approximately that of glass). If outside the
		 * sphere overlap, an index of 1.000 is returned (approximately
		 * that of air).
		 */
		inline
		double
		nuValue
			( engabra::g3::Vector const & rVec
			) const
		{
			double nu{ 1.000 }; // default values is near that of air

			// define a double convex lens as intrsection of two spheres
			std::function<bool(engabra::g3::Vector const &)> const inLens
				{ [] (engabra::g3::Vector const & rVec)
					{
						using namespace engabra::g3;
						constexpr double r1Sq{ sq(11.) };
						static Vector const c1{ -10., 0., 0. };
						constexpr double r2Sq{ sq(21.) };
						static Vector const c2{  20., 0., 0. };
						bool const in1{ (magSq(rVec - c1) < r1Sq) };
						bool const in2{ (magSq(rVec - c2) < r2Sq) };
						return (in1 && in2);
					}
				};

			if (inLens(rVec))
			{
				nu = 1.500; // inside the "lens" set IoR to be glass-like
			}
			return nu;
		}

	}; // Media
	// [DoxyExampleMedia]
}


/*! \brief Provides a complete example of AeroPlygiant use.
 *
 * \include demoAeroPlygiant.cpp
 */
int
main
	()
{
	// Define an active volume that is of interest to the situation at hand.
	std::shared_ptr<aply::env::ActiveVolume> const ptVolume
		{ std::make_shared<aply::env::ActiveBox>
			( engabra::g3::Vector{ -5., -10., -10. }
			, engabra::g3::Vector{  5.,  10.,  10. }
			)
		};
	// Specity this volume as a clipping region in the base class
	// (ref aply::env::IndexVolume::qualifiedNuValue() method).
	Media const media(ptVolume);

	// Specify initial conditions (tangent direction and first point on path)
	using engabra::g3::Vector;
	Vector const tanBeg{ direction(Vector{ 1., .2, .3 }) }; // tangent dir
	Vector const locBeg{ -5., 0., 0. }; // first point on ray
	aply::ray::Start const start{ aply::ray::Start::from(tanBeg, locBeg) };

	// Configure propagation step size and specify path save interval
	constexpr double propStepDist{ 1./1024. };
	constexpr double saveStepDist{ 1./16. };

	// An approximate end point can be used in ray::Path ctor to
	// estimate and allocate path storage space. This is useful if
	// ray path is nominally follows a smooth "kind of straight" curve.
	// Alternatively, can use ray::Path::reserve() function to
	// explicitly allocate a specific amount of space (in which case
	// the ctor default argument approxEndLoc=null<Vector>() is used.
	Vector const approxEndLoc{ 10., 0., 0. };
	aply::ray::Path path(start, saveStepDist, approxEndLoc);

	// Create propagor engine and request it to trace path
	// Same "prop" instance can be used for many paths.
	aply::ray::Propagator const prop{ &media, propStepDist };
	prop.tracePath(&path);

	// Access individual nodes in traced path.
	for (aply::ray::Node const & node : path.theNodes)
	{
		// Here, just display brief summary of individual node information
		std::cout << node.infoBrief() << std::endl;
	}

	return 0;
}

