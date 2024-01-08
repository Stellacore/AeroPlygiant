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

#ifndef aply_examp_DiffEqSystem_INCL_
#define aply_examp_DiffEqSystem_INCL_

/*! \file
 *
 * \brief Example classes for differential equation systems.
 *
 */


#include "mathDiffEqSystem.hpp"

#include <Engabra>


namespace aply
{
namespace examp
{

/*! \brief Differential equation related examples
 *
 */
namespace diffeq
{

	/*! \brief System of equations associated with uniform acceleration.
	 *
	 * This system of equations models a uniform acceleration such as
	 * what is experienced by an object falling in the gravity field
	 * near the surface of the Earth (where gravity acceleration is
	 * nominally uniform and constant).
	 *
	 * The relevant differential equation is:
	 * \arg y'' = g   (for constant gravity acceleration value 'g')
	 *
	 * The associated simultaneous equation system structure is
	 * described in functor operator()().
	 *
	 * Initial values are described in method initValues().
	 */
	struct UniformAccel : public aply::math::DiffEqSystem
	{
		//! Initial Time: Associated to initial conditions
		double const theTime0{ engabra::g3::null<double>() };

		//! Initial Height: from which object is dropped (e.g. y at t0)
		double const theHeight0{ engabra::g3::null<double>() };

		//! Initial Speed: with which object is moving (e.g. y' at t0)
		double const theSpeed0{ engabra::g3::null<double>() };

		//! Nominal value of local gravity (near 45-deg lattitude 1k' elevation)
		static constexpr double theAccel{ -9.805 };

		//! No-op construction
		explicit
		UniformAccel
			( double const & t0  //!< Time [sec] at which drop starts
			, double const & h0  //!< Height [m] at which drop starts
			, double const & v0  //!< Speed [m/sec] with which drop starts
			)
			: DiffEqSystem()
			, theTime0{ t0 }
			, theHeight0{ h0 }
			, theSpeed0{ v0 }
		{ }

		//! Default (no-op) dtor
		virtual
		~UniformAccel
			() = default;

		/*! \brief Derivative equation system function values.
		 *
		 * The relevant parameter and functions are:
		 * \arg t:  evolution parameter - e.g. time
		 * \arg y0 = y:  position function (expect: y = g*t*t/2 + y0*t + v0)
		 * \arg y1 = y0' = y': velocity function (expect: y = g*t + y0)
		 * \arg y2 = y1' = y'': acceleration (expect: y = g)
		 *
		 * This functor, operator(), gets input argument values:
		 * \arg t = xyValues.first  (e.g. time)
		 * \arg y0 = xyValues.second[0] (e.g. position)
		 * \arg y1 = xyValues.second[1] (e.g. velocity)
		 *
		 * This functor, operator(), produces output function values:
		 * \arg y0' = y1
		 * \arg y1' = g + (0.*y1 + 0.*y0)  // example doesn't depend on y0, y1
		 */
		// [DoxyExampleDiffEqSystemOp()]
		virtual
		std::vector<double>
		operator()
			( std::pair<double, std::vector<double> > const & xyValues
				//!< values represent [(x), (y0, y1, ...)]
			) const
		{
			double const & xValue = xyValues.first;
			std::vector<double> const & yFuncs = xyValues.second;
			double const & y0Func = yFuncs[0];
			double const & y1Func = yFuncs[1];

			double const y0Prime{ y1Func };
			double const y1Prime{ theAccel };
			return std::vector<double>
				{ y0Prime
				, y1Prime
				};
		}
		// [DoxyExampleDiffEqSystemOp()]

		/*! \brief Init values: evoluation parameter, position, and velocity.
		 *
		 * Return value is:
		 * \code
		 * return
		 *	{ theTime0    // time zero
		 * 	, { theHeight // position (y0) at theTime0
		 * 	  , 0.        // velocity (y1) at theTime0
		 * 	  }
		 *	};
		 * \endcode
		 */
		// [DoxyExampleDiffEqSystemInitVal]
		virtual
		std::pair<double, std::vector<double> >
		initValues
			() const
		{
			double const t0{ theTime0 };
			double const y0_t0{ theHeight0 };
			double const y1_t0{ theSpeed0 };
			return { t0, std::vector<double>{ y0_t0, y1_t0 } };
		}
		// [DoxyExampleDiffEqSystemInitVal]

		//! Expected acceleration at time tau (from known analytical solution)
		double
		expAccelerationAt
			( double const & tau
			) const
		{
			return theAccel;
		}

		//! Expected velocity at time tau (from known analytical solution)
		double
		expVelocityAt
			( double const & tau
			) const
		{
			double const dTau{ tau - theTime0 };
			return
				( theAccel * dTau
				+ theSpeed0
				);
		}

		//! Expected position at time tau (from known analytical solution)
		double
		expPositionAt
			( double const & tau
			) const
		{
			double const dTau{ tau - theTime0 };
			return
				( .5 * theAccel * dTau * dTau
				+ theSpeed0 * dTau
				+ theHeight0
				);
		}


	}; // UniformAccel

} // [diffeq]

} // [examp]
} // [aply]

#endif // aply_examp_DiffEqSystem_INCL_

