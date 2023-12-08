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

#include <Engabra>



namespace tst
{
	using namespace engabra::g3;

	/*! \brief Thick slab of constant index of refraction
	 *
	 * Classic "thick plate" refraction model
	 *
	 */
	struct Slab : public env::IndexVolume
	{
		Vector const theNormDir{ null<Vector>() };
		double const theBegDot{ null<double>() };
		double const theEndDot{ null<double>() };
		double const theNuPrev{ null<double>() };
		double const theNuCurr{ null<double>() };
		double const theNuNext{ null<double>() };

		//! Value constructor
		inline
		explicit
		Slab
			( Vector const & normDir
			, double const & begDot
			, double const & endDot
			, double const & nuPrev = 1.
			, double const & nuCurr = 1.5
			, double const & nuNext = 1.
			)
			: IndexVolume{}
			, theNormDir{ direction(normDir) }
			, theBegDot{ begDot }
			, theEndDot{ endDot }
			, theNuPrev{ nuPrev }
			, theNuCurr{ nuCurr }
			, theNuNext{ nuNext }
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
			double nu{ null<double>() };
			double const valDot{ (rVec * theNormDir).theSca[0] };
			if (valDot < theBegDot)
			{
				nu = theNuPrev;
			}
			else
			if ((theBegDot < valDot) && (valDot < theEndDot))
			{
				nu = theNuCurr;
			}
			else
			if (theEndDot < valDot)
			{
				nu = theNuNext;
			}
			return nu;
		}

	}; // Slab

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


	//! Atmospheric model : nu = alpha*exp(-beta*radius)
	struct AtmModel : public env::IndexVolume
	{
		//! Classic exponential decay model
		inline
		static
		double
		expValue
			( double const & amp
			, double const & scale
			, double const & radius
			)
		{
			return amp * std::exp(-scale * radius);
		}

		//! Exponential decay constant
		inline
		static
		double
		beta
			( double const & v1
			, double const & v2
			, double const & r1
			, double const & r2
			)
		{
			return { log(v1 / v2) / (r2 - r1) };
		}

		//! Exponential decay magnitude
		inline
		static
		double
		alpha
			( double const & v1
			, double const & v2
			, double const & r1
			, double const & r2
			)
		{
			double const bNeg{ - beta(v1, v2, r1, r2) };
			return
				{ (v1 + v2) / (std::exp(bNeg * r1) + std::exp(bNeg * r2)) };
		}


		double const theRadGround{ null<double>() };
		double const theRadSpace{ null<double>() };
		double const theAlpha{ null<double>() };
		double const theBeta{ null<double>() };

		//! Construct an invalid instance
		inline
		AtmModel
			()
			: IndexVolume{}
			, theAlpha{ null<double>() }
			, theBeta{ null<double>() }
		{
		}

		//! Construct model to match environment constants
		inline
		AtmModel
			( env::Planet const & planet
			)
			: theRadGround{ planet.theRadGround }
			, theRadSpace{ planet.theRadSpace }
			, theAlpha
				{ alpha
					( planet.theNuGround
					, planet.theNuSpace
					, theRadGround
					, theRadSpace
					)
				}
			, theBeta
				{ beta
					( planet.theNuGround
					, planet.theNuSpace
					, theRadGround
					, theRadSpace
					)
				}
		{ }

		//! Thickness of atmosphere
		inline
		double
		thickness
			() const
		{
			return theRadSpace - theRadGround;
		}

		//! Index of refraction value at vector location rVec
		virtual
		inline
		double
		nuValue
			( Vector const & rVec
			) const
		{
			double const rMag{ magnitude(rVec) };
			return expValue(theAlpha, theBeta, rMag);
		}

		//! Gradient (approximation) of Index of refraction at location rVec
		virtual
		inline
		Vector
		nuGradient
			( Vector const & rVec
			) const
		{
			double const mag{ magnitude(rVec) };
			double const del
				{ mag * std::sqrt(std::numeric_limits<double>::epsilon()) };
			double const dNu1
				{ nuValue(rVec + del * e1) - nuValue(rVec - del * e1) };
			double const dNu2
				{ nuValue(rVec + del * e2) - nuValue(rVec - del * e2) };
			double const dNu3
				{ nuValue(rVec + del * e3) - nuValue(rVec - del * e3) };
			double const scale{ 1. / (2. * del) };
			return { scale * Vector{ dNu1, dNu2, dNu3 } };
		}

		//! Sampling of atm.nuValue() values
		inline
		std::vector<double>
		nuProfile
			( double const & delta
			, double const & rBeg = env::sEarth.theRadGround
			, double const & rEnd = env::sEarth.theRadSpace
			) const
		{
			std::vector<double> nus;
			double const dSteps{ (rEnd - rBeg) / delta };
			std::size_t const numSamps{ static_cast<std::size_t>(dSteps) +2u };
			nus.reserve(numSamps);
			for (double rad{rBeg} ; rad < rEnd ; rad += delta)
			{
				// evaluate in arbitrary direction (here e3)
				double const nu{ nuValue(rad * e3) };
				nus.emplace_back(nu);
			}
			return nus;
		}

	}; // AtmModel


} // [tst]

#endif // Refraction_tstModels_INCL_

