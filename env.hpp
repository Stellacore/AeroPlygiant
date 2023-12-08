//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//

#ifndef Refraction_env_INCL_
#define Refraction_env_INCL_

/*! \file
 *
 * \brief Environment configuration parameters (related to Refraction)
 *
 */


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

} // [env]

#endif // Refraction_env_INCL_

