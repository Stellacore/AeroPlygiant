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
 * \brief Unit test for class Hessian
 *
 */


#include "tst.hpp"
#include "math.hpp"

#include <Engabra>

#include <cmath>
#include <limits>
#include <sstream>
#include <vector>


namespace
{
	using namespace engabra::g3;

	//! \brief A scalar field with spatial variation for use in testing.
	struct ScalarField
	{
		//! Arbitrary function for testing
		inline
		double
		operator()
			( Vector const & loc
			) const
		{
			double const & xx = loc[0];
			double const & yy = loc[1];
			double const & zz = loc[2];
			return
				(  3. * xx * xx * xx
				+  5. * yy * yy * yy
				+  7. * zz * zz * zz
				+ 11. * xx * yy
				+ 13. * yy * zz
				+ 17. * zz * xx
				+ 23.
				);
		}

		double dfdx( Vector const & loc ) const
			{ return 3.*3.*loc[0]*loc[0] + 11.*loc[1] + 17.*loc[2]; };
		double dfdy( Vector const & loc ) const
			{ return 5.*3.*loc[1]*loc[1] + 11.*loc[0] + 13.*loc[2]; };
		double dfdz( Vector const & loc ) const
			{ return 7.*3.*loc[2]*loc[2] + 13.*loc[1] + 17.*loc[0]; };

		double dfdxdx( Vector const & loc ) const { return 3.*3.*2.*loc[0]; };
		double dfdxdy( Vector const & loc ) const { return 11.; };
		double dfdxdz( Vector const & loc ) const { return 17.; };

		double dfdydx( Vector const & loc ) const { return 11.; };
		double dfdydy( Vector const & loc ) const { return 5.*3.*2.*loc[1]; };
		double dfdydz( Vector const & loc ) const { return 13.; };

		double dfdzdx( Vector const & loc ) const { return 17.; };
		double dfdzdy( Vector const & loc ) const { return 13.; };
		double dfdzdz( Vector const & loc ) const { return 7.*3.*2.*loc[2]; };

		inline
		aply::math::Matrix
		hessian
			( Vector const & loc
			) const
		{
			return aply::math::Matrix
				{ { dfdxdx(loc), dfdxdy(loc), dfdxdz(loc) }
				, { dfdydx(loc), dfdydy(loc), dfdydz(loc) }
				, { dfdzdx(loc), dfdzdy(loc), dfdzdz(loc) }
				};
		}

	}; // ScalarField


	//! Check second order derivative components
	void
	test0
		( std::ostringstream & oss
		)
	{
		// [DoxyExample00]
		
		// Location at which to evaluate function
		Vector const loc{ 10., 20., 30. };

		// Function/Functor that generates scalar value at location
		// (and here, provides exact analytica values for test comparison)
		ScalarField const func;

		// Compute numeric estimate of gradient
		Vector const gradF{ aply::math::gradientOf(func, loc) };

		// Numerically estimated Hessian
		aply::math::Matrix const gotHess{ aply::math::hessianOf(func, loc) };

		// [DoxyExample00]

		// expected values
		aply::math::Matrix const expHess{ func.hessian(loc) };

		double const tol
			{ std::sqrt(std::numeric_limits<double>::epsilon()) };
		if (! aply::math::nearlyEquals(gotHess, expHess, tol))
		{
			aply::math::Matrix const difHess{ gotHess - expHess };
			oss << "Failure of hessian test" << std::endl;
			oss << "exp:\n" << expHess << std::endl;
			oss << "got:\n" << gotHess << std::endl;
			oss << "dif:\n" << difHess << std::endl;
		}
	}
}


/*! \brief Unit test for Hessian component computation
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);

	return tst::finish(oss);
}

