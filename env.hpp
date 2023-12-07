//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//

#ifndef Refraction_env_INCL_
#define Refraction_env_INCL_

//! \file Environment configuration parameters (related to Refraction)


#include <Engabra>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>


/*! \brief Functions and data associated with ray propagation environment
 */
namespace env
{
	using namespace engabra::g3;

	//! Data parameters for Planet atmosphere
	struct Planet
	{
		double const theNuGround{ null<double>() };
		double const theNuSpace{ null<double>() };
		double const theRadGround{ null<double>() };
		double const theRadSpace{ null<double>() };

		//! Descriptive information about this instance
		inline
		std::string
		infoString
			( std::string const & title = {}
			) const
		{
			std::ostringstream oss;
			if (! title.empty())
			{
				oss << title << '\n';
			}
			oss << " theNuGround: " << io::fixed(theNuGround);
			oss << '\n';
			oss << "  theNuSpace: " << io::fixed(theNuSpace);
			oss << '\n';
			oss << "theRadGround: " << io::fixed(theRadGround);
			oss << '\n';
			oss << " theRadSpace: " << io::fixed(theRadSpace);
			return oss.str();
		}

	}; // Planet

	//! Parameters for Earth atmosphere
	static Planet const sEarth
		{
		  1.000273 // nuGround // (for earth, actual is 1.000273 at STP)
		, 1.000    // nuSpace
		, 6370.e3  // radiusGround
		, 6470.e3  // radiusSpace
		};

	//! Interface specification for refractive media volume
	struct IndexVolume
	{
		//! Index of refraction value at vector location rVec
		inline
		virtual
		double
		nuValue
			( Vector const & rVec
			) const = 0;

		/*! \brief Gradient (approximate) for Index of refraction at rVec.
		 *
		 * Default implementation estimates gradient numerically.
		 */
		inline
		virtual
		Vector
		nuGradient
			( Vector const & rVec
				//!< Location at which to estimate gradient
			, double const & stepSize
				//!< Use half this value as difference estimating gradient
			) const
		{
			double const del{ .5 * stepSize };
			double const scl{ 1. / stepSize };
			return Vector
				{ scl * (nuValue(rVec + del*e1) - nuValue(rVec - del*e1))
				, scl * (nuValue(rVec + del*e2) - nuValue(rVec - del*e2))
				, scl * (nuValue(rVec + del*e3) - nuValue(rVec - del*e3))
				};
		}

	}; // IndexVolume


	//! Atmospheric model : nu = alpha*exp(-beta*radius)
	struct AtmModel : public IndexVolume
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
			: theAlpha{ null<double>() }
			, theBeta{ null<double>() }
		{
		}

		//! Construct model to match environment constants
		inline
		AtmModel
			( Planet const & planet
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
			, double const & rBeg = sEarth.theRadGround
			, double const & rEnd = sEarth.theRadSpace
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

} // [env]

#endif // Refraction_env_INCL_

