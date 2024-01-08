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


#include "envAtmosphereParameters.hpp"
#include "mathDiffEqSolve.hpp"
#include "rayRefraction.hpp"

#include <sstream>


namespace
{

	//! \brief Refraction solver system of equations.
	class RefractionSystem : public aply::math::DiffEqSystem
	{
	public:

		aply::ray::Refraction const & theRefraction;

		RefractionSystem
			( aply::ray::Refraction const & refraction
			)
			: theRefraction(refraction)
		{
		}

		virtual
		std::vector<double>
		operator()
			( std::pair<double, std::vector<double> > const & in
			) const
		{
			std::vector<double> out;

			double const pointRadius(in.first);

			double const & refractiveInvariant
				(theRefraction.theRefractiveInvariant);
			double const pointElevation
				(pointRadius - theRefraction.theRadiusEarth);
			double const refraction
				(theRefraction.theAtmosphere.indexOfRefraction(pointElevation));
			using engabra::g3::sq;
			double const radicand
				(sq(pointRadius*refraction) - sq(refractiveInvariant));

			out.push_back
				(refractiveInvariant / pointRadius / std::sqrt(radicand));

			return out;
		}

		virtual
		std::pair<double, std::vector<double> >
		initValues
			() const
		{
			return theRefraction.theInitValues;
		}
	};

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
	, theInitValues(std::make_pair(startHeight, std::vector<double>(1, 0.0)))
{
}

bool
Refraction :: isValid() const
{
	return engabra::g3::isValid(theRadiusEarth);
}

double
Refraction :: angleAt
	( double const & radius
	) const
{
	math::DiffEqSolve solver(50.0); // reasonable step
	theInitValues = solver.solutionFor(radius, RefractionSystem(*this));
	return theInitValues.second[0];
}

double
Refraction :: displacementAt
	( double const & radius
	) const
{
	return radius * angleAt(radius);
}

} // [ray]
} // [aply]

