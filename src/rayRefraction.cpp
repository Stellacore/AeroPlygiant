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
		std::pair<double, std::vector<double> > const theInitValues
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
			, std::pair<double, std::vector<double> > const & initValues
			, aply::env::Atmosphere const & atmosphere
			, double const & radEarth
			)
			: DiffEqSystem{}
			, theRefConst{ refConst }
			, theInitValues{ initValues }
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
			return theInitValues;
		}

	}; // RefractGyer

} // [anon]


namespace aply
{
namespace ray
{

Refraction :: Refraction()
	: theRadiusEarth{ engabra::g3::null<double>() }
	, theAtmosphere()
	, theRefractiveInvariant(0.0)
	, theInitValues()
{
}

Refraction :: Refraction
	( double const & radiusEarth
	, double const & startHeight
	, double const & startAngle
	)
	: theRadiusEarth(radiusEarth)
	, theAtmosphere(env::Atmosphere::COESA1976())
	, theRefractiveInvariant
		( startHeight * theAtmosphere.indexOfRefraction(startHeight-radiusEarth)
		* std::sin(startAngle))
	, theInitValues
		{ std::make_pair(startHeight, std::vector<double>{ 0. }) }
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
	math::DiffEqSolve solver(50.0); // reasonable step
	RefractGyer const refractionSystem
		( theRefractiveInvariant
		, theInitValues
		, theAtmosphere
		, theRadiusEarth
		);
	std::pair<double, std::vector<double> > const endValues
		{ solver.solutionFor(radius, refractionSystem) };
	return endValues.second[0];
}

double
Refraction :: displacementAt
	( double const & radius
	) const
{
	return radius * thetaAngleAt(radius);
}

} // [ray]
} // [aply]

