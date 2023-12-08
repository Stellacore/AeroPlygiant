//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


#ifndef Refraction_tstModels_INCL_
#define Refraction_tstModels_INCL_


/*! \file
 * 
 * \brief Index of refraction volume models for testing and demonstration.
 *
 */


#include "env.hpp"
// #include "ray.hpp"
// #include "save.hpp"

#include <Engabra>

// #include <iostream>
// #include <sstream>



namespace tst
{
	using namespace engabra::g3;

	/*! \brief Simple example of a spherical shape with constant index.
	 */
	struct Sphere : public env::IndexVolume
	{
		Vector const theCenter{ null<Vector>() };
		double const theRadius{ null<double>() };
		double const theNuCenter{ null<double>() };
		double const theNuEdge{ null<double>() };

		//! Construct a sphere in space
		inline
		explicit
		Sphere
			( Vector const & center
			, double const & radius
			, double const & nuCenter = 1.5
			, double const & nuEdge = 1.
			)
			: IndexVolume{}
			, theCenter{ center }
			, theRadius{ radius }
			, theNuCenter{ nuCenter }
			, theNuEdge{ nuEdge }
		{
		}

		//! Index of refraction value at vector location rVec
		inline
		virtual
		double
		nuValue
			( Vector const & rVec
			) const
		{
			double nu{ theNuEdge };
			// assume linear gradient
			// nu(r) = nuCenter - frac*(nuEdge-nuCenter)
			// nu(r) = nuCenter - (r/radius)*(nuEdge-nuCenter)
			// nu(r) = nuCenter - r*(1./radius)*(nuEdge-nuCenter)
			// grad(nu) = -(1./radius)*nuEdge-nuCenter;
			double const dist{ magnitude(rVec - theCenter) };
			double const frac{ dist / theRadius };
			if (frac < 1.)
			{
				nu = frac*(theNuEdge - theNuCenter) + theNuCenter;
			}
			return nu;
		}

		//! Override Gradient (approximation) with analytical expression
		inline
		virtual
		Vector
		nuGradient
			( Vector const & rVec
			, double const & //
			) const
		{
			Vector grad{ zero<Vector>() };
			Vector const delta{ rVec - theCenter };
			double const dist{ magnitude(delta) };
			if (dist < theRadius)
			{
				Vector const gDir{ direction(delta) };
				double const gMag{ (-1./theRadius)*(theNuEdge-theNuCenter) };
				grad = gMag * gDir;
			}
			return grad;
		}

	}; // IndexVolume


	//! \brief Thick slab of constant index of refraction
	struct Slab : public env::IndexVolume
	{
		double const theBegX{ null<double>() };
		double const theEndX{ null<double>() };
		double const theNuInside{ null<double>() };
		double const theNuOutside{ null<double>() };

		//! Value constructor
		inline
		explicit
		Slab
			( double const & begX
			, double const & endX
			, double const & nuInside = 1.5
			, double const & nuOutside = 1.
			)
			: IndexVolume{}
			, theBegX{ begX }
			, theEndX{ endX }
			, theNuInside{ nuInside }
			, theNuOutside{ nuOutside }
		{
		}


		//! \brief Index of refraction value at vector location rVec.
		inline
		virtual
		double
		nuValue
			( Vector const & rVec
			) const
		{
			double nu{ theNuOutside };
			double const & valX = rVec[0];
			if ((theBegX < valX) && (valX < theEndX))
			{
				nu = theNuInside;
			}
			return nu;
		}

	}; // Slab

} // [tst]

#endif // Refraction_tstModels_INCL_

