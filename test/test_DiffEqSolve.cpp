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
 *
 * \brief Unit test for class math::DiffEqSolve
 *
 */


#include "mathDiffEqSolve.hpp"
#include "mathDiffEqSystem.hpp"

#include "../example/diffeqSystem.hpp"
#include "tst.hpp"

#include <Engabra>

#include <sstream>


namespace
{
#if 0
	/*! \brief Example system of equations associated with a falling stone.
	 *
	 * This system of equations models a object falling in gravity field
	 * near the surface of the Earth (where gravity field can be 
	 * considered uniform and constant).
	 *
	 * The relevant differential equation is:
	 * \arg y'' = g   (for constant gravity acceleration value 'g')
	 */
	struct RockDrop : public aply::math::DiffEqSystem
	{
		//! Initial Time: Associated to initial conditions
		double const theTime0{ engabra::g3::null<double>() };

		//! Initial Height: from which object is dropped (e.g. y at t0)
		double const theHeight0{ engabra::g3::null<double>() };

		//! Initial Speed: with which object is moving (e.g. y' at t0)
		double const theSpeed0{ engabra::g3::null<double>() };

		//! Nominal value of local gravity (near 45-deg lattitude 1k' elevation)
		static constexpr double theAccelG{ -9.805 };

		//! No-op construction
		explicit
		RockDrop
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
		~RockDrop
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
			double const y1Prime{ theAccelG };
			return std::vector<double>
				{ y0Prime
				, y1Prime
				};
		}

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

	}; // RockDrop
#endif


	//! Check if numeric solution agrees with known analytical solution.
	void
	test0
		( std::ostringstream & oss
		)
	{
		// [DoxyExample00]

		// Construct "rock drop" equations of motion system

		// initial conditions
		constexpr double dummyOffset{ 100. };
		constexpr double t0{  0. + dummyOffset };
		constexpr double h0{ 10. - dummyOffset };
		constexpr double v0{  0. };

		// end point of integration for this test
		constexpr double t1{  2. + t0 };

		// specify equation system
		aply::examp::diffeq::UniformAccel const uniAccelEquations(t0, h0, v0);

		// configure integrator
		constexpr double stepSize{ .001 };
		aply::math::DiffEqSolve solver(stepSize);

		// request solution at time t1
		std::pair<double, std::vector<double> > const soln
			{ solver.solutionFor(t1, uniAccelEquations)
			};

		// interpret returned solution data
		double const & got_t1 = soln.first;  // should match t1
		std::vector<double> const & yVals = soln.second;
		double const & gotPos = yVals[0];
		double const & gotVel = yVals[1];

		// known analytical solution
		double const expAcc{ uniAccelEquations.expAccelerationAt(t1) };
		double const expVel{ uniAccelEquations.expVelocityAt(t1) };
		double const expPos{ uniAccelEquations.expPositionAt(t1) };


		// utilize engabra capabilities for compare
		using engabra::g3::nearlyEquals;
		using engabra::g3::io::fixed;
		using engabra::g3::io::enote;

		if (! nearlyEquals(got_t1, t1))
		{
			oss << "Failure of RockDrop end time (t1) test\n";
			oss << "exp_t1: " << fixed(    t1) << '\n';
			oss << "got_t1: " << fixed(got_t1) << '\n';
		}

		constexpr double tol{ 1.e-12 };
		if (! nearlyEquals(gotPos, expPos, tol))
		{
			double const difPos{ gotPos - expPos };
			oss << "Failure of RockDrop position test\n";
			oss << "expPos: " << fixed(expPos) << '\n';
			oss << "gotPos: " << fixed(gotPos) << '\n';
			oss << "difPos: " << enote(difPos) << '\n';
		}
		if (! nearlyEquals(gotVel, expVel, tol))
		{
			double const difVel{ gotVel - expVel };
			oss << "Failure of RockDrop position test\n";
			oss << "expVel: " << fixed(expVel) << '\n';
			oss << "gotVel: " << fixed(gotVel) << '\n';
			oss << "difVel: " << enote(difVel) << '\n';
		}

		if (! oss.str().empty())
		{
			oss << '\n';
			oss << "     exp @t1: " << fixed(t1) << '\n';
			oss << "     exp Pos: " << fixed(expPos) << '\n';
			oss << "     exp Vel: " << fixed(expVel) << '\n';
			oss << "     exp Acc: " << fixed(expAcc) << '\n';
			oss << '\n';
			oss << " soln got_t1: " << fixed(got_t1) << '\n';
			oss << "soln yVals #: " << yVals.size() << '\n';
			oss << "soln yVal[0]: " << fixed(yVals[0]) << '\n';
			oss << "soln yVal[1]: " << fixed(yVals[1]) << '\n';
			oss << '\n';
		}
		// [DoxyExample00]
	}

} // [anon]


/*! \brief Unit test for math::DiffEqSolve
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss); // Falling object equation of motion test

	return tst::finish(oss);
}

