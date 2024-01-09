//
// Original Code:
//  Copyright (c) 2007 Stellacore Corporation.
//  Donated to AeroPlygiant Open Source project 2024.
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
\brief Definitions for ray::Refraction
*/


#include "mathDiffEqSolve.hpp"
#include "rayRefraction.hpp"

#include <sstream>


namespace
{

	/*! \brief Implementation of Gyer paper Eqn[13] integration.
 	*/
	struct RefractGyer : public aply::math::DiffEqSystem
	{
		//! \brief Refraction constant (invariant along ray).
		double const theRefConst{ engabra::g3::null<double>() };

		/*! \brief Initial conditions comprising....
		 *
		 *	- .first: starting height (relative to Earth \b center)
		 *	- .second: vector of size one
		 *		- : constant of integration for Theta_c
		 *
		 */
		std::pair<double, std::vector<double> > const theInitRadTheta
			{ std::make_pair
				( engabra::g3::null<double>()
				, std::vector<double>{}
				)
			};

		//! \brief Atmosphere model in location of interest.
		aply::env::Atmosphere const theAtmosphere{ };

		//! \brief Radius of Earth in vicinity of location of interest.
		double const theRadEarth{ engabra::g3::null<double>() };


		//! \brief Value ctor.
		explicit
		RefractGyer
			( double const & refConst
			, std::pair<double, std::vector<double> > const & initRadTheta
			, aply::env::Atmosphere const & atmosphere
			, double const & radEarth
			)
			: DiffEqSystem{}
			, theRefConst{ refConst }
			, theInitRadTheta{ initRadTheta }
			, theAtmosphere{ atmosphere }
			, theRadEarth{ radEarth }
		{ }

		//! \brief No-op dtor.
		virtual
		~RefractGyer
			() = default;


		/*! \brief Single ODE from Gyer Eqn[12].
		 *
		 * Implements integration of Eqn(12) in Gyer's paper
		 *   https://www.asprs.org/wp-content/uploads/
		 *     pers/1996journal/mar/1996_mar_301-310.pdf
		 *
		 * I.e.
		 * \arg Theta_c = integral{ k / (r*sqrt(n*n*r*r - k*k)) * dr };
		 */
		virtual
		std::vector<double>
		operator()
			( std::pair<double, std::vector<double> > const & in
			) const
		{
			std::vector<double> derivFuncs(1u);

			// Current evaluation functions
			double const currRad(in.first);
			// std::vector<double> const & currRnFuncs = in.second;

			// height relative to Earth radius
			double const elev{ currRad - theRadEarth };
			double const currIoR{ theAtmosphere.indexOfRefraction(elev) };
			using engabra::g3::sq;
			double const radicand{ sq(currRad*currIoR) - sq(theRefConst) };
			double const denom{ currRad * std::sqrt(radicand) };

			// Derivative function values
			double const r0Prime{ theRefConst / denom }; // integrand

			// Derivative function values
			return std::vector<double>
				{ r0Prime
				};
		}

		//! \brief Start height and inital 'Theta_c' value (generally 0.)
		virtual
		std::pair<double, std::vector<double> >
		initValues
			() const
		{
			return theInitRadTheta;
		}

	}; // RefractGyer

	/*! \brief Info on net ray deviation as observed from sensor station.
 	 */
	struct NetRayInfo
	{
		//! Distance from \b center of Earth at which ray starts.
		double const theBegRadius{ engabra::g3::null<double>() };

		//! Viewing angle from Nadir direction (0. is straight down).
		double const theBegLookAngle{ engabra::g3::null<double>() };


		//! \brief Deviation (refracted w.r.t. ideal straight line) at Sensor.
		double
		refractionDeviation
			( double const endRadius
				//!< Distance from \b center of Earth at which ray terminates.
			, double const endTheta
				//!< Angle subtended by ray path \b from \b Earth \b center.
			) const
		{
			using namespace engabra::g3;

			// [DoxyExample01]

			// relative to local Nadir topocentric frame
			static Vector const upDir{ e3 };
			static Vector const downDir{ -upDir };

			// Cartesian locations in local polar frame
			Vector const locBeg{ theBegRadius * upDir };
			Vector const locEnd
				{ endRadius * std::sin(endTheta)
				, 0.
				, endRadius * std::cos(endTheta)
				};
			Vector const locDel{ locEnd - locBeg };

			// compute angular deviation at the sensor (in local Nadir frame)
			BiVector const deviationAngle3D{ (logG2(downDir * locDel)).theBiv };
			double const deviationAngle{ magnitude(deviationAngle3D) };
			double const deviationMag{ theBegLookAngle - deviationAngle };

			// [DoxyExample01]
			return deviationMag;
		}

	}; // NetRayInfo

} // [anon]


namespace aply
{
namespace ray
{

Refraction :: Refraction()
	: theRadiusEarth{ engabra::g3::null<double>() }
	, theAtmosphere()
	, theRefractiveInvariant(0.0)
	, theInitRadTheta()
{
}

Refraction :: Refraction
	( double const & lookAngle
	, double const & radiusSensor
	, double const & radiusEarth
	, env::Atmosphere const & atmosphere
	)
	: theStartLookAngle{ lookAngle }
	, theStartRadius{ radiusSensor }
	, theRadiusEarth{ radiusEarth }
	, theAtmosphere{ atmosphere }
	, theRefractiveInvariant
		{ theStartRadius
		* theAtmosphere.indexOfRefraction(theStartRadius-theRadiusEarth)
		* std::sin(lookAngle)
		}
	, theInitRadTheta
		{ std::make_pair(theStartRadius, std::vector<double>{ theTheta0 }) }
{
}

bool
Refraction :: isValid() const
{
	return engabra::g3::isValid(theRadiusEarth);
}

double
Refraction :: thetaAngleAt
	( double const & radius
	) const
{
	// For aerial sensing work, an integration step size of 50 [m] seems
	// to be a good value. E.g. 10x larger or 10x smaller still produces
	// the same ray deviation angle from 9k[m] at pi/4 look dir.
	constexpr double stepSize{ 50. }; // a resonable step size

	math::DiffEqSolve solver(stepSize);
	RefractGyer const refractionSystem
		( theRefractiveInvariant
		, theInitRadTheta
		, theAtmosphere
		, theRadiusEarth
		);
	// Return initial value like structure includes:
	//	- endValues.first : radius (should match method input argument)
	//	- endValues.second: size of 1u
	//		- [0] : Theta_c angle (ray path polar angle from center of Earth)
	std::pair<double, std::vector<double> > const endValues
		{ solver.solutionFor(radius, refractionSystem) };
	return endValues.second[0];
}

double
Refraction :: angularDeviationFromStart
	( double const radiusEnd
	, double const thetaEnd
	) const
{
	double deviation{ engabra::g3::null<double>() };
	if (isValid())
	{
		NetRayInfo const netRayInfo{ theStartRadius, theStartLookAngle };
		deviation = netRayInfo.refractionDeviation(radiusEnd, thetaEnd);
	}
	return deviation;
}

} // [ray]
} // [aply]

