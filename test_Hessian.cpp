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

#include <Engabra>

#include <cmath>
#include <limits>
#include <sstream>
#include <vector>


namespace math
{
	//! Numeric estimate of gradient for a scalar function
	template <typename Func>
	inline
	engabra::g3::Vector
	gradientOf
		( Func const & func
		, engabra::g3::Vector const & loc
		, double const & relStepSize
			= std::sqrt(std::numeric_limits<double>::epsilon())
		)
	{
		using namespace engabra::g3;
		double const fval{ func(loc) };
		double const stepSize{ relStepSize * fval };
		double const del{ .5 * stepSize };
		double const scl{ 1. / stepSize };
		return Vector
			{ scl * ( func(loc + del*e1) - func(loc - del*e1) )
			, scl * ( func(loc + del*e2) - func(loc - del*e2) )
			, scl * ( func(loc + del*e3) - func(loc - del*e3) )
			};
	}

	using Row = std::vector<double>;
	using Matrix = std::vector<Row>;

	template <typename Type>
	inline
	Type
	null
		()
	{
		using engabra::g3::null;
		return Type
			{ { null<double>(), null<double>(), null<double>() }
			, { null<double>(), null<double>(), null<double>() }
			, { null<double>(), null<double>(), null<double>() }
			};
	}

	//! Numeric estimate of Hessian for a scalar function
	template <typename Func>
	inline
	Matrix
	hessianOf
		( Func const & func
		, engabra::g3::Vector const & loc
		, double const & relStepSize
			= std::sqrt(std::numeric_limits<double>::epsilon())
		)
	{
		using namespace engabra::g3;
		Matrix hess{ null<Matrix>() };
		double const delta{ .1 };
		Vector const d1{ delta * e1 };
		Vector const d2{ delta * e2 };
		Vector const d3{ delta * e3 };
		double const scale{ 1. / sq(delta) };
		{ // dxdy = dydx
			double const aa{ func(loc-d1-d2) };
			double const bb{ func(loc+d1-d2) };
			double const cc{ func(loc-d1+d2) };
			double const dd{ func(loc+d1+d2) };
			double const dxdy{ scale*(.25 * (dd - bb - cc + aa)) }; // dydx
			hess[0][1] = dxdy;
			hess[1][0] = dxdy;
		}
		{ // dydz = dzdy
			double const aa{ func(loc-d2-d3) };
			double const bb{ func(loc+d2-d3) };
			double const cc{ func(loc-d2+d3) };
			double const dd{ func(loc+d2+d3) };
			double const dydz{ scale*(.25 * (dd - bb - cc + aa)) }; // dzdy
			hess[1][2] = dydz;
			hess[2][1] = dydz;
		}
		{ // dzdx = dxdz
			double const aa{ func(loc-d3-d1) };
			double const bb{ func(loc+d3-d1) };
			double const cc{ func(loc-d3+d1) };
			double const dd{ func(loc+d3+d1) };
			double const dzdx{ scale*(.25 * (dd - bb - cc + aa)) }; // dxdz
			hess[2][0] = dzdx;
			hess[0][2] = dzdx;
		}
		{ // dxdx
			double const aa{ func(loc-d1) };
			double const bb{ func(loc) };
			double const dd{ func(loc+d1) };
			double const dxdx{ scale*((dd - 2.*bb + aa)) };
			hess[0][0] = dxdx;
		}
		{ // dydy
			double const aa{ func(loc-d2) };
			double const bb{ func(loc) };
			double const dd{ func(loc+d2) };
			double const dydy{ scale*((dd - 2.*bb + aa)) };
			hess[1][1] = dydy;
		}
		{ // dzdz
			double const aa{ func(loc-d3) };
			double const bb{ func(loc) };
			double const dd{ func(loc+d3) };
			double const dzdz{ scale*((dd - 2.*bb + aa)) };
			hess[2][2] = dzdz;
		}

		return hess;
	}

} // [math]


namespace
{
	inline
	bool
	isValid
		( math::Matrix const & matrix
		)
	{
		bool valid{ false };
		if (3u == matrix.size())
		{
			valid = 
				(  (3u == matrix[0].size())
				&& (3u == matrix[1].size())
				&& (3u == matrix[2].size())
				);
		}
		return valid;
	}

	inline
	bool
	nearlyEquals
		( math::Matrix const & got
		, math::Matrix const & exp
		, double const & tol = std::numeric_limits<double>::epsilon()
		)
	{
		bool same{ false };
		if (isValid(got) && isValid(exp))
		{
			same = true; // until/unless otherwise
			for (std::size_t ii{0u} ; ii < 3u ; ++ii)
			{
				for (std::size_t jj{0u} ; jj < 3u ; ++jj)
				{
					using engabra::g3::nearlyEquals;
					same &= nearlyEquals(got[ii][jj], exp[ii][jj], tol);
					if (! same)
					{
						break;
					}
				}
			}
		}
		return same;
	}

	inline
	std::ostream &
	operator<<
		( std::ostream & ostrm
		, math::Matrix const & matrix
		)
	{
		std::size_t const numDigPre{ 5u };
		std::size_t const numDigPost{ 9u };
		using engabra::g3::io::fixed;
		if (3u == matrix.size())
		{
			if (3u << matrix[0].size())
			{
				ostrm
					<< ' ' << fixed(matrix[0][0], numDigPre, numDigPost)
					<< ' ' << fixed(matrix[0][1], numDigPre, numDigPost)
					<< ' ' << fixed(matrix[0][2], numDigPre, numDigPost)
					;
			}
			if (3u << matrix[1].size())
			{
				ostrm << '\n';
				ostrm
					<< ' ' << fixed(matrix[1][0], numDigPre, numDigPost)
					<< ' ' << fixed(matrix[1][1], numDigPre, numDigPost)
					<< ' ' << fixed(matrix[1][2], numDigPre, numDigPost)
					;
			}
			if (3u << matrix[2].size())
			{
				ostrm << '\n';
				ostrm
					<< ' ' << fixed(matrix[2][0], numDigPre, numDigPost)
					<< ' ' << fixed(matrix[2][1], numDigPre, numDigPost)
					<< ' ' << fixed(matrix[2][2], numDigPre, numDigPost)
					;
			}
		}
		return ostrm;
	}

	inline
	math::Matrix
	operator-
		( math::Matrix const & matA
		, math::Matrix const & matB
		)
	{
		math::Matrix result{ math::null<math::Matrix>() };
		for (std::size_t ii{0u} ; ii < 3u ; ++ii)
		{
			for (std::size_t jj{0u} ; jj < 3u ; ++jj)
			{
				result[ii][jj] = matA[ii][jj] - matB[ii][jj];
			}
		}
		return result;
	}

}


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
		math::Matrix
		hessian
			( Vector const & loc
			) const
		{
			return math::Matrix
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
		Vector const gradF{ math::gradientOf(func, loc) };

		// Numerically estimated Hessian
		math::Matrix const gotHess{ math::hessianOf(func, loc) };

		// [DoxyExample00]

		// expected values
		math::Matrix const expHess{ func.hessian(loc) };

		double const tol
			{ std::sqrt(std::numeric_limits<double>::epsilon()) };
		if (! nearlyEquals(gotHess, expHess, tol))
		{
			math::Matrix const difHess{ gotHess - expHess };
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

