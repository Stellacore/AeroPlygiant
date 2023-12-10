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

#ifndef aply_env_IndexVolume_INCL_
#define aply_env_IndexVolume_INCL_

/*! \file
 *
 * \brief Environment configuration parameters (related to Refraction)
 *
 */


#include "envActiveVolume.hpp"

#include <Engabra>

#include <memory>


/*! \brief Functions and data associated with ray propagation environment
 */
namespace env
{
	using namespace engabra::g3;

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

#endif // aply_env_IndexVolume_INCL_

