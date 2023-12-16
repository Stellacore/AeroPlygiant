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
 * \brief Unit test for algebraic equation solution (used in Refraction.lyx)
 *
 */


#include "tst.hpp"

#include <Engabra>

#include <iostream>


namespace
{
	using namespace engabra::g3;

	//! Evaluate left-hand-side (a*x-x*b)
	inline
	MultiVector
	lhs
		( MultiVector const & mv_a
		, MultiVector const & mv_b
		, MultiVector const & mv_x
		)
	{
		return mv_a * mv_x + mv_x * mv_b;
	}

	//! Evaluate full equation (a*x-x*b - D)
	inline
	MultiVector
	equation
		( MultiVector const & mv_a
		, MultiVector const & mv_b
		, MultiVector const & mv_x
		, MultiVector const & mv_D
		)
	{
		return (lhs(mv_a, mv_b, mv_x) - mv_D);
	}
}


/*! \brief Numeric cross check on solution to sylvesters equation.
 */
int
main
	()
{
	std::ostringstream oss;

	using namespace engabra::g3;

	// create a bivector equation of the form encountered in Refraction DiffEq
	Vector const gVec{ 2., 3., 5. };
	Vector const tVec{ 11., 13., 17. };

	MultiVector const mv_a{ tVec + gVec };
	MultiVector const mv_b{ tVec - gVec };
	MultiVector const mv_D{ BiVector{ 19., 23., 27. } };

	MultiVector const expEqn{ zero<MultiVector>() };
	double const tol
		{ magnitude(mv_D) * std::numeric_limits<double>::epsilon() };

	/*
	std::cout << "mv_a: " << mv_a << '\n';
	std::cout << "mv_b: " << mv_b << '\n';
	std::cout << "mv_D: " << mv_D << '\n';
	*/

	// check with MultiVector arithmetic

	double const aSq{ magSq(mv_a) };
	double const bSq{ magSq(mv_b) };
	MultiVector const mv_aInv{ (1./aSq) * mv_a };
	MultiVector const mv_aDvr{ dirverse(mv_a) };
	MultiVector const mv_bDvr{ dirverse(mv_b) };

	// // General solution form (for vector coefficients and biv RHS
	// MultiVector const coef{ mv_a + mv_aInv*mv_b*mv_bDvr + mv_b + mv_bDvr };

	//
	// using a^-1 in derivation
	//

	{
	MultiVector const coef{ mv_a - mv_aInv*mv_b*mv_b };
	MultiVector const fact{ mv_D - mv_aInv*mv_D*mv_b };
	MultiVector const cInv{ inverse(coef.theVec) }; // okay here
	MultiVector const soln{ cInv * fact };

	MultiVector const gotEqn{ equation(mv_a, mv_b, soln, mv_D) };
	tst::checkGotExp(oss, gotEqn, expEqn, "soln(1a)", tol);
	}

	{
	MultiVector const coef{ mv_aInv * (mv_a * mv_a - mv_b * mv_b) };
	MultiVector const fact{ mv_aInv * (mv_a * mv_D - mv_D * mv_b) };

	MultiVector const cInv{ inverse(coef.theVec) }; // okay here
	MultiVector const soln{ cInv * fact };

	MultiVector const gotEqn{ equation(mv_a, mv_b, soln, mv_D) };
	tst::checkGotExp(oss, gotEqn, expEqn, "soln(1b)", tol);
	}

	{
	MultiVector const fact{ (mv_a * mv_D - mv_D * mv_b) };
	MultiVector const soln{ (1./(aSq-bSq)) * fact };

	MultiVector const gotEqn{ equation(mv_a, mv_b, soln, mv_D) };
	tst::checkGotExp(oss, gotEqn, expEqn, "soln(1c)", tol);
	}

	{
	MultiVector const soln{ (1./(aSq-bSq)) * (mv_a*mv_D - mv_D*mv_b) };

	MultiVector const gotEqn{ equation(mv_a, mv_b, soln, mv_D) };
	tst::checkGotExp(oss, gotEqn, expEqn, "soln(1d)", tol);
	}

	{
	double const scl{ 1./(aSq-bSq) };
	MultiVector const sum{ mv_a + mv_b };
	MultiVector const dif{ mv_a - mv_b };
	MultiVector const vec{ scl * (sum * mv_D).theVec };
	MultiVector const tri{ scl * (dif * mv_D).theTri };
	MultiVector const soln{ vec + tri };

	MultiVector const gotEqn{ equation(mv_a, mv_b, soln, mv_D) };
	tst::checkGotExp(oss, gotEqn, expEqn, "soln(1e)", tol);
	}

	// using b^-1 in derivation

	{
	double const scl{ 1./(bSq-aSq) };
	MultiVector const top{ mv_D * mv_b - mv_a * mv_D };
	MultiVector const soln{ scl * top };

	MultiVector const gotEqn{ equation(mv_a, mv_b, soln, mv_D) };
	tst::checkGotExp(oss, gotEqn, expEqn, "soln(2a)", tol);
	}

	return tst::finish(oss);
}

