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

#ifndef aply_INCL_
#define aply_INCL_

/*! \file
 *
 * \brief Environment configuration parameters (related to Refraction)
 *
 */


#include <Engabra>

#include <cmath>
#include <vector>


namespace aply
{

/*! \brief Mathematical utility functions
 *
 */
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

	//! \brief A 1D row of a grid data structure
	using Row = std::vector<double>;
	//! \brief A 2D gridded data structure
	using Grid = std::vector<Row>;

	//! Matrix as a grid of data values.
	using Matrix = Grid;

	//! Null type for Matrix i.e. use: null<Matrix>()
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

	//! True if Matrix has the correct size
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
		// and first element is not null
		if (valid)
		{
			valid &= engabra::g3::isValid(matrix[0][0]);
		}
		return valid;
	}

	//! True if got and exp are element by element engabra::g3::nearlyEquals()
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

	//! \brief Subtraction of two matrices. Result = A-B
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

} // [math]


} // [aply]

//! global scope overloads
namespace
{
	//! Put aply::math::Matrix to standard out
	inline
	std::ostream &
	operator<<
		( std::ostream & ostrm
		, aply::math::Matrix const & matrix
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

	//! Publish matrix operator
	using aply::math::operator-;

} // [anon]

#endif // aply_INCL_


