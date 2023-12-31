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
	constexpr
	double
	mBarForPascal
		( double const & pPascal
		)
	{
		return .01 * pPascal;
	}


	//! \brief Kelvin for degrees Celsius.
	inline
	constexpr
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
	struct CylindricalAir : public aply::env::IndexVolume
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
		CylindricalAir
			( aply::geom::Cylinder const & tube
			, double const & tempOnAxisK
			, double const & tempOnEdgeK
			)
			: IndexVolume{}
			, theTube{ tube }
			, theNuInterval{ nuInterval(tempOnAxisK, tempOnEdgeK) }
		{ }

		//! Index of refraction associated with radial gradient along cylinder
		inline
		double
		nuValue
			( engabra::g3::Vector const & rLoc
			) const
		{
			double nu{ engabra::g3::null<double>() };
			double const lenFrac{ theTube.fractionAlongAxis(rLoc) };
			if ((! (lenFrac < 0.)) && (lenFrac < 1.))
			{
				nu = theNuInterval.max(); // default to STP air
				double const radFrac{ theTube.fractionFromAxis(rLoc) };
				if (radFrac < 1.)
				{
					nu = theNuInterval.valueAtFrac(radFrac);
				}
			}
			return nu;
		}

	}; // CylindricalAir

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
	constexpr double hotRadius{ 10. }; // [m]
	constexpr double endPad{ 1. };
	constexpr double length{ 250. + endPad }; // [m]

	constexpr double tempOnAxisK{ units::kelvinForC(35.) };
	constexpr double tempOnEdgeK{ units::kelvinForC(25.) };

	static Vector const axisDir{ e2 };
	static Vector const offsetDir{ e1 };
	static Vector const elevDir{ e3 };
	static Vector const axisBeg{ -endPad * axisDir };
	static Vector const axisEnd{ length * axisDir };

	static Vector const staLoc{                - 5.*offsetDir + 1.5*elevDir };
	static Vector const tgtLoc{ length*axisDir + 0.*offsetDir + 1.5*elevDir };

	static double const obsDist{ magnitude(tgtLoc - staLoc) };
	static Vector const obsDir{ direction(tgtLoc - staLoc) };

	static Vector const approxEndLoc{ length * obsDir };

	// tracing configuration
	constexpr double propStepDist{ .001 }; // Propagation step size[m]
	constexpr double saveStepDist{ 10. }; // Data save step size[m]

	using namespace aply;

	// Configure a cylinder with axis centered on 'road' center
	road::CylindricalAir const media
		{ geom::Cylinder
			( axisBeg    // axis starting at origin
			, axisDir    // axis along the positive 'x' direction
			, obsDist    // length [m]
			, hotRadius  // radius [m]
			)
		, tempOnAxisK
		, tempOnEdgeK
		};

	// setup ray start parallel to cylinder 'off to the side' of the 'road'
	ray::Start const start{ ray::Start::from(obsDir, staLoc) };

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

