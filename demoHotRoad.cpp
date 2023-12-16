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


namespace units
{
	//! \brief MillBars for value expressed in Pascal.
	inline
	double
	mBarForPascal
		( double const & pPascal
		)
	{
		return .01 * pPascal;
	}


	//! \brief Kelvin for degrees Celsius.
	inline
	double
	kelvinForC
		( double const & degC
		)
	{
		return 273.15 + degC;
	}

} // [units]

namespace air
{
	//
	// Standard conditions:
	// -- https://en.wikipedia.org/wiki/Standard_temperature_and_pressure
	//

	//! Standard temperature [K]
	constexpr double sStdTemperature{ 293.15 };

	//! Standard pressure [Pa]
	constexpr double sStdPressure{ 101325. };

	//! Standard relative humidity [fraction]
	constexpr double sStdRelHumidity{ 0.00 };

	/*! \brief Index of Refraction for given temperature and pressure.
	 *
	 * Formula from Gyer 1996 (ref Paper.bib) eqn (13).
	 *
	 * Site:
	 *   https://refractiveindex.info/
	 * Python script for evaluating Ciddor equation:
	 *   https://github.com/polyanskiy/refractiveindex.info-scripts/
	 *   blob/master/scripts/Ciddor%201996%20-%20air.py
	 *
	 * Sample values from: https://emtoolbox.nist.gov/Wavelength/Ciddor.asp
	 *
	 * # Temp(C)  Pres(kPa)   IoR
	 *
	 *	-20    100.000    1.000310769  *uncertain
	 *	  0    100.000    1.000287830
	 *	 20    100.000    1.000267817
	 *	 40    100.000    1.000249811
	 *
	 *	-20     80.000    1.000248567  *uncertain
	 *	  0     80.000    1.000230213
	 *	 20     80.000    1.000214152
	 *	 40     80.000    1.000199587
	 *
	 *	-20     60.000    1.000186388
	 *	  0     60.000    1.000172610
	 *	 20     60.000    1.000160495
	 *	 40     60.000    1.000149367
	 */
	inline
	double
	nuForTP
		( double const & airTempK
			//!< Temperature in Kelvin
		, double const & airPresPa
			//!< Pressure in Pascals (Newton per m^2)
		)
	{
		double const & mBarPres{ .01 * airPresPa };
		// formula from Gyer1996
		double const refractivity{ .000078831 * (mBarPres / airTempK) };
		double const nu{ 1. + refractivity };
		return nu;
	}

} // [air]


//! Utilities for supporting HotRoad demo
namespace road
{
	/*! \brief Cylindrical volume refactive index varying by radius.
	 *
	 * Intended to represent the changing index of refraction such as
	 * due to hot air accumulating above and around a long straight road.
	 *
	 */
	struct AirVolume : public aply::env::IndexVolume
	{
		//! Cylindrical tube of (linearly) varying air IoR
		aply::geom::Cylinder const theTube;
		//! Index of refraction gap from axis to outside radial edge
		aply::geom::Interval const theNuInterval;

		/*! \brief An index of refraction gradient in radial direction.
		 *
		 * Index of refraction is estimated based on provided air
		 * temperatures at center (on axis) and edge of the cylinder.
		 *
		 * IoR formula extracted from Gyer 1996 PE&RS article. (Ref
		 * Papers.bib).
		 */
		inline
		static
		aply::geom::Interval
		nuInterval
			( double const & airTempOnAxis
			, double const & airTempAtEdge
			)
		{
			double const nuAxis
				{ air::nuForTP(airTempOnAxis, air::sStdPressure) };
			double const nuEdge
				{ air::nuForTP(airTempAtEdge, air::sStdPressure) };
			return aply::geom::Interval(nuAxis, nuEdge);
		}

		//! Construct this shape and alignment
		inline
		explicit
		AirVolume
			( aply::geom::Cylinder const & tube
			)
			: IndexVolume{}
			, theTube{ tube }
			, theNuInterval
				{ nuInterval
					( units::kelvinForC(35.)
					, units::kelvinForC(25.)
					)
				}
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

	static Vector const staLoc{           - 5.*e2 };
	static Vector const tgtLoc{ length*e1 + 5.*e2 };

	static Vector const obsDir{ direction(tgtLoc - staLoc) };
	static Vector const approxEndLoc{ length * obsDir };

	// tracing configuration
	constexpr double propStepDist{ .001 }; // Propagation step size[m]
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
		{ ray::Start::from(obsDir, staLoc + 1.e-6*e1) };

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

