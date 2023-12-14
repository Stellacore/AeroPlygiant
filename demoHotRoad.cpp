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
 * \brief Demonstration: ray path adacent to low index air bubble.
 *
 */


#include "geom.hpp"
#include "env.hpp"
#include "ray.hpp"

#include <iostream>
#include <utility>


//! Utilities for supporting HotRoad demo
namespace road
{
	//! \brief TODO
	struct AirVolume : public aply::env::IndexVolume
	{
		//! Cylindrical tube of (linearly) varying air IoR
		aply::geom::Cylinder const theTube;
		//! Index of refraction gap from axis to outside radial edge
		aply::geom::Interval const theNuInterval;

		//! Create (fake) index of refraction gradient in radial direction
		inline
		static
		aply::geom::Interval
		nuInterval
			()
		{
			double const nuGround{ aply::env::sEarth.theNuGround };
			double const refractivity{ nuGround - 1. }; // refractivity (IoR-1.)
			double const nuAxis{ 1. + .5*refractivity };
			return aply::geom::Interval(nuAxis, nuGround);
		}

		//! Construct this shape and alignment
		inline
		explicit
		AirVolume
			( aply::geom::Cylinder const & tube
			)
			: IndexVolume{}
			, theTube{ tube }
			, theNuInterval{ nuInterval() }
		{ }

		//! Index of refraction associated with radial gradient along cylinder
		inline
		double
		nuValue
			( engabra::g3::Vector const & rLoc
			) const
		{
			double nu{ theNuInterval.max() }; // default to STP air
			double const lenFrac{ theTube.fractionAlongAxis(rLoc) };
			if ((! (lenFrac < 0.)) && (lenFrac < 1.))
			{
				double const radFrac{ theTube.fractionFromAxis(rLoc) };
				if (radFrac < 1.)
				{
					nu = theNuInterval.valueAtFrac(radFrac);
				}
			}
			return nu;
		}

	}; // AirVolume

} // [road]


/*! \brief Simulate survey sighting along the edge of a hot roadway.
 *
 * Hot air above roadway is simulated with a half-cylinder IndexVolume
 * shape aligned with the road. Refractivity at the center is 1/2 that
 * at the edge (which is ambient air STP index).
 *
 * The propagated path starts aligned with the cylinder, at 1/2 radius
 * form the axis. The ray is propagated forward to the end of the cylinder
 * and the associated path curvature results are reported.
 *
 */
int
main
	( int argc
	, char * argv[]
	)
{
	using namespace engabra::g3;

	// scene configuration
	static Vector const axisDir{ e1 };
	constexpr double hotRadius{ 10. }; // [m]
	constexpr double length{ 1000. }; // [m]
	static Vector const staDir{ -e2 }; // direction to station from origin
	static Vector const obsDir{ axisDir }; // observation direction
	static Vector const approxEndLoc{ length * obsDir };

	// tracing configuration
//	constexpr double propStepDist{ .001 }; // Propagation step size[m]
//	constexpr double saveStepDist{ .001 }; // Data save step size[m]
constexpr double propStepDist{  10. }; // Propagation step size[m]
constexpr double saveStepDist{ 100. }; // Data save step size[m]

	using namespace aply;

	// Configure a cylinder with axis centered on 'road' center
	road::AirVolume const media
		{ geom::Cylinder
			( zero<Vector>()  // axis starting at origin
			, axisDir         // axis along the positive 'x' direction
			, length          // length [m]
			, hotRadius       // radius [m]
			)
		};

// TODO
// NOTE: the initial index of refraction is starting at "air" value
//       (somewhere inside Propagator maybe?)
//       should take starting NU value from environment.

	// setup ray start parallel to cylinder 'off to the side' of the 'road'
	ray::Start const start
		{ ray::Start::from(obsDir, .5*hotRadius*staDir + 1.e-6*e1) };

	// construct propagator
	ray::Propagator const prop{ &media, propStepDist };

	// trace path
	ray::Path path(start, saveStepDist, approxEndLoc);
	prop.tracePath(&path);

	// report results
	for (ray::Node const & node : path.theNodes)
	{
		std::cout << node.infoBrief() << '\n';
	}
	std::cout << ray::PathView(&path).infoCurvature() << std::endl;

}

