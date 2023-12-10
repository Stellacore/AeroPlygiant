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
	using namespace aply;


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
			, double const & nuPrev
			, double const & nuCurr
			, double const & nuNext
			, std::shared_ptr<env::ActiveVolume>
				const & ptVolume = env::sPtAllSpace
			)
			: IndexVolume(ptVolume)
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
			, std::shared_ptr<env::ActiveVolume>
				const & ptVolume = env::sPtAllSpace
			)
			: IndexVolume(ptVolume)
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


	//! \brief An exponential decay function that matches boundary values.
	struct ExpDecay
	{
		//! Decay constant - matching v0(r0) and v1(r1)
		inline
		static
		double
		beta
			( double const & v0
			, double const & v1
			, double const & r0
			, double const & r1
			)
		{
			return { log(v0 / v1) / (r1 - r0) };
		}

		//! Amplitude factor - matching v0(r0) and v1(r1)
		inline
		static
		double
		lnAlpha
			( double const & v0
			, double const & v1
			, double const & r0
			, double const & r1
			)
		{
			double const frac{ 1. / (r1 - r0) };
			return (r1 * frac * std::log(v0) - r0 * frac * std::log(v1));
		}

		//! Amplitude factor - matching v0(r0) and v1(r1)
		inline
		static
		double
		alpha
			( double const & v0
			, double const & v1
			, double const & r0
			, double const & r1
			)
		{
			return std::exp(lnAlpha(v0, v1, r0, r1));
		}

		double const theAlpha{ null<double>() }; //!< Amplitude factor
		double const theBeta{ null<double>() }; //!< Decay constant (mag)

		//! Construct default null instance.
		ExpDecay
			()
			: theAlpha{ null<double>() }
			, theBeta{ null<double>() }
		{ }

		/*! \brief An exponential decay function matching boundary values.
		 *
		 * The resulting function (via the value() method) provides a
		 * value that decays exponentially as:
		 * \arg value = alpha * exp(-beta * rad)
		 *
		 * with:
		 * \arg v0 = value(r0)
		 * \arg v1 = value(r1)
		 */
		explicit
		ExpDecay
			( double const & v0
			, double const & v1
			, double const & r0
			, double const & r1
			)
			: theAlpha{ alpha(v0, v1, r0, r1) }
			, theBeta{ beta(v0, v1, r0, r1) }
		{ }

		//! Classic exponential decay model
		inline
		double
		operator()
			( double const & radius
			) const
		{
			return theAlpha * std::exp(-theBeta * radius);
		}

		//! Descriptive information about this instance.
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
			oss << "theAlpha: " << io::fixed(theAlpha)
				<< " "
				<< "theBeta: " << io::fixed(theBeta)
				;
			return oss.str();
		}


	}; // ExpDecay


	//! Atmospheric model : nu = alpha*exp(-beta*radius)
	struct AtmModel : public env::IndexVolume
	{
		std::pair<double, double> const the_v0r0{}; //!< 1st boundary loc/val
		std::pair<double, double> const the_v1r1{}; //!< 2nd boundary loc/val
		ExpDecay const theNuFunc{};

		//! Construct an invalid instance
		inline
		AtmModel
			()
			: IndexVolume()
			, theNuFunc{}
		{ }

		//! Construct model to match environment constants
		inline
		AtmModel
			( env::Planet const & planet
			, std::shared_ptr<env::ActiveVolume>
				const & ptVolume = env::sPtAllSpace
			)
			: IndexVolume(ptVolume)
			, the_v0r0{ planet.theNuGround, planet.theRadGround }
			, the_v1r1{ planet.theNuSpace, planet.theRadSpace }
			, theNuFunc
				( the_v0r0.first
				, the_v1r1.first
				, the_v0r0.second
				, the_v1r1.second
				)
		{ }

		//! Thickness of atmosphere
		inline
		double
		thickness
			() const
		{
			return (the_v1r1.second - the_v0r0.second);
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
			return theNuFunc(rMag);
		}

		//! Descriptive information about this instance.
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

			oss << "v0: " << io::fixed(the_v0r0.first)
				<< " "
				<< "at r0: " << io::fixed(the_v0r0.second)
				;
			oss << '\n';
			oss << "v1: " << io::fixed(the_v1r1.first)
				<< " "
				<< "at r1: " << io::fixed(the_v1r1.second)
				;
			oss << '\n';
			oss << theNuFunc.infoString();
			return oss.str();
		}

	}; // AtmModel


} // [tst]

#endif // Refraction_tstModels_INCL_

