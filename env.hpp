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

#ifndef Refraction_env_INCL_
#define Refraction_env_INCL_

/*! \file
 *
 * \brief Environment configuration parameters (related to Refraction)
 *
 */


#include <Engabra>

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>


/*! \brief Functions and data associated with ray propagation environment
 */
namespace env
{
	using namespace engabra::g3;

	/*! \brief Data parameters for spherical planet's atmosphere.
	 *
	 * Container for data describing useful for creating a simple
	 * spherically symmetric index of refraction field (e.g. for
	 * describing Earth atmospher models)
	 *
	 * The 'nu' naming refer to Index of Refraction (IoR) values.
	 */
	struct Planet
	{

		//! Index of refraction ('nu') at ground level.
		double const theNuGround{ null<double>() };
		//! Index of refraction ('nu') at edge of space
		double const theNuSpace{ null<double>() };

		//! Radius from planet center with IoR == theNuGround
		double const theRadGround{ null<double>() };
		//! Radius from planet center with IoR == theNuSpace
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

	/*! \brief Specify volume of space through which rays should be propagated.
	 */
	struct ActiveVolume
	{
		std::string theName{};

		//! Construct a named instance
		explicit
		ActiveVolume
			( std::string const & name = "ActiveVolume"
			)
			: theName{ name }
		{ }

		//! Overload to define shape of volume (true: inside, false: outside)
		inline
		virtual
		bool
		contains
			( Vector const & rVec
			) const
		{
std::cout << "ActiveVolume::contains" << std::endl;
			return true;
		}

	}; // ActiveVolume

	//! An active volume w/o limits
	static std::shared_ptr<ActiveVolume> const sPtAllSpace
		{ std::make_shared<ActiveVolume>("sAllSpace") };

	/*! \brief Interface specification for refractive media volume
	 *
	 * Represent the ray trace environment as a volume of refractive
	 * index values. The value can vary arbitrarily in order to
	 * simulate continuous changes in index of refraction or discrete
	 * shapes.
	 * 
	 * \note the boundaries (of interest) are indicated by returning
	 * a null value for refraction index from method nuValue().
	 *
	 */
	struct IndexVolume
	{
		/*! \brief Region in which ray propagation should be performed.
		 *
		 * By default, the index volume (IoR field) is active everywhere.
		 * E.g. ray propagation will never hit an edge
		 */
		std::shared_ptr<ActiveVolume> const thePtVolume{};

		//! \brief Construct media IoR volume (clipped by ActiveVolume)
		inline
		explicit
		IndexVolume
			( std::shared_ptr<ActiveVolume> const & ptVolume = sPtAllSpace
			)
			: thePtVolume{ ptVolume }
		{ }

		//! \brief Index of Refraction value, or null if outside thePtVolume.
		inline
		double
		qualifiedNuValue
			( Vector const & rVec
			) const
		{
			double nu{ null<double>() }; // default to stop condition
			if (thePtVolume->contains(rVec))
			{
				nu = nuValue(rVec);
			}
std::cout
	<< "Volume: " << thePtVolume->theName
	<< "  testing rVec: " << rVec
	<< "  nu: " << nu
	<< std::endl;

			return nu;
		}

		/*! \brief Index of refraction value at vector location rVec.
		 *
	 	 * Note: return nuValue = null<double>() to indicate the edges
		 * of the volume (boundaries at which ray tracing operations
		 * should stop.
		 */
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

