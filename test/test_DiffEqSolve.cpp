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

	//! Check if numeric solution agrees with known analytical solution.
	void
	test0
		( std::ostringstream & oss
		)
	{
		// Test with uniform acceleration equation system

		// initial conditions
		constexpr double dummyOffset{ 100. };
		constexpr double t0{  0. + dummyOffset };
		constexpr double h0{ 10. - dummyOffset };
		constexpr double v0{  0. };

		// end point of integration for this test
		constexpr double t1{  2. + t0 };

		// [DoxyExample00]

		// specify equation system (inherits from aply::math::DiffEqSystem)
		aply::examp::diffeq::UniformAccel const uniAccelEquations(t0, h0, v0);

		// configure integrator
		constexpr double stepSize{ .001 };
		aply::math::DiffEqSolve solver(stepSize);

		// request solution at time t1
		// solution includes: { time, {funcVal, func'Val, func''Val, ...} }
		std::pair<double, std::vector<double> > const soln
			{ solver.solutionFor(t1, uniAccelEquations)
			};

		// interpret returned solution data
		// for uniform acceleration equation system,
		// ... second integrated function is position (funcVal)
		// ... first integrated function is velocity (funcVal')
		// ... original func (acceleration) from UniformAccel (not returned)
		double const & got_t1 = soln.first;  // should match t1
		std::vector<double> const & yVals = soln.second;
		double const & gotPos = yVals[0];
		double const & gotVel = yVals[1];

		// [DoxyExample00]

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
			oss << "Failure of UniformAccel end time (t1) test\n";
			oss << "exp_t1: " << fixed(    t1) << '\n';
			oss << "got_t1: " << fixed(got_t1) << '\n';
		}

		constexpr double tol{ 1.e-12 };
		if (! nearlyEquals(gotPos, expPos, tol))
		{
			double const difPos{ gotPos - expPos };
			oss << "Failure of UniformAccel position test\n";
			oss << "expPos: " << fixed(expPos) << '\n';
			oss << "gotPos: " << fixed(gotPos) << '\n';
			oss << "difPos: " << enote(difPos) << '\n';
		}
		if (! nearlyEquals(gotVel, expVel, tol))
		{
			double const difVel{ gotVel - expVel };
			oss << "Failure of UniformAccel position test\n";
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

