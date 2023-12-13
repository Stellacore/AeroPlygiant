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


#include "env.hpp"
#include "ray.hpp"

#include <iostream>
#include <utility>


//! \brief Basic geometry primatives.
namespace geom
{
	using namespace engabra::g3;

	/*! \brief Use to values to define a distance scale (origin and unit value)
	 *
	 * Perhaps best explained by example:
	 * \arg fracAtValue(1., pair<>(2., 3.)) = -1 // extrapolation
	 * \arg fracAtValue(2., pair<>(2., 3.)) =  0 // Begin interval Included
	 * \arg fracAtValue(3., pair<>(2., 3.)) =  1 // End interval EXcluded
	 * \arg fracAtValue(4., pair<>(2., 3.)) =  2 // extrapolation
	 *
	 * And the inverses:
	 * \arg valueAtFrac(.75, pair<>(2., 3.)) == 2.75
	 *
	 */
	struct Interval
	{
		//! Define the half open interval [min,max)
		std::pair<double, double> const theMinMax{};
		double const theSpan{ null<double>() };
		double const theScale{ null<double>() };

		explicit
		Interval
			( double const & begValue
			, double const & endValue
			)
			: theMinMax{ begValue, endValue }
			, theSpan{ theMinMax.second - theMinMax.first }
			, theScale{ 1. / theSpan }
		{ }

		//! The origin of the interval
		inline
		double
		min
			() const
		{
			return theMinMax.first;
		}

		//! The origin of the interval
		inline
		double
		max
			() const
		{
			return theMinMax.second;
		}

		/*! \brief Inter(extra)polated fraction of way into interval.
		 */
		inline
		double
		fracAtValue
			( double const & value
			) const
		{
			return theScale * (value - min());
		}

		//! \brief Value associated with fraction between end points.
		inline
		double
		valueAtFrac
			( double const & frac
			) const
		{
			return (frac * theSpan + min());
		}

	}; // Interval

	//! \brief A geometric cylinder shape of finite length
	struct Cylinder
	{
		Vector const theAxisBeg;
		Vector const theAxisDir;
		double const theLength;
		double const theRadius;
		Interval const theGapLength;
		Interval const theGapRad;

		//! Value construction
		inline
		explicit
		Cylinder
			( Vector const & axisBeg
			, Vector const & axisDir
			, double const & length
			, double const & radius
			)
			: theAxisBeg{ axisBeg }
			, theAxisDir{ direction(axisDir) }
			, theLength{ length }
			, theRadius{ radius }
			, theGapLength{ 0., theLength }
			, theGapRad{ 0., theRadius }
		{ }

		//! Distance orthogonal from body axis to loc.
		inline
		double
		distanceFromAxis
			( Vector const & someLoc
			) const
		{
			double const dist
				{ ((someLoc - theAxisBeg)*theAxisDir).theSca[0] };
			return dist;
		}

		//! Distance orthogonal from body axis to loc.
		inline
		double
		fractionFromAxis
			( Vector const & someLoc
			) const
		{
			return theGapRad.fracAtValue(distanceFromAxis(someLoc));
		}

		//! Distance parallel along body axis to loc.
		inline
		double
		distanceAlongAxis
			( Vector const & someLoc
			) const
		{
			double const dist
				{ (someLoc*theAxisDir).theSca[0] };
			return dist;
		}

		inline
		double
		fractionAlongAxis
			( Vector const & someLoc
			) const
		{
			return theGapLength.fracAtValue(distanceAlongAxis(someLoc));
		}

	}; // Cylinder

} // [geom]

//! Utilities for supporting HotRoad demo
namespace road
{
	//! \brief TODO
	struct AirVolume : public aply::env::IndexVolume
	{
		geom::Cylinder const theTube; //!< cylindrical tube of reduced air IoR
		geom::Interval const theGapNu;

		//! Construct this shape and alignment
		inline
		explicit
		AirVolume
			( geom::Cylinder const & tube
			)
			: IndexVolume{}
			, theTube{ tube }
			, theGapNu
				{ .5 * aply::env::sEarth.theNuGround  // IoR along axis
				,      aply::env::sEarth.theNuGround  // IoR along outer edge
				}
		{ }

		//! Index of refraction associated with radial gradient along cylinder
		inline
		double
		nuValue
			( engabra::g3::Vector const & rLoc
			) const
		{
			double nu{ theGapNu.max() }; // default to STP air
			double const lenFrac{ theTube.fractionAlongAxis(rLoc) };
			if ((! (lenFrac < 0.)) && (lenFrac < 1.))
			{
				double const radFrac{ theTube.fractionFromAxis(rLoc) };
				if (radFrac < 1.)
				{
					nu = theGapNu.valueAtFrac(radFrac);
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
	static Vector const dirToRight{ -e2 };
	constexpr double hotRadius{ 10. }; // [m]

	// tracing configuration
	constexpr double propStepDist{ .001 }; // Propagation step size[m]
	constexpr double saveStepDist{ .001 }; // Data save step size[m]

	// Configure a cylinder with axis centered on 'road' center
	road::AirVolume const media
		{ geom::Cylinder
			( zero<Vector>()  // axis starting at origin
			, axisDir         // axis along the positive 'x' direction
			, 1000.           // length [m]
			, hotRadius       // radius [m]
			)
		};

	// setup ray start parallel to cylinder 'off to the side' of the 'road'
	using namespace aply;
	ray::Start const start
		{ ray::Start::from(axisDir, .5*hotRadius*dirToRight) };

	// construct propagator
	ray::Propagator const prop{ &media, propStepDist };

	// trace path
	ray::Path path(start, saveStepDist);
	prop.tracePath(&path);

	// report results
	std::cout << ray::PathView(&path).infoCurvature() << std::endl;

}

