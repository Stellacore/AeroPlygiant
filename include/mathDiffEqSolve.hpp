//
// Original Code:
//  Copyright (c) 2007 Stellacore Corporation.
//  Donated to AeroPlygiant Open Source project 2024.
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

#ifndef aply_math_DiffEqSolve_INCL_
#define aply_math_DiffEqSolve_INCL_

/*! \file
\brief Declarations for math::DiffEqSolve
*/


#include "mathDiffEqSystem.hpp"

#include <string>
#include <vector>


namespace aply
{
namespace math
{

/*! \brief This class solves ordinary differential equations numerically.

The ODE equation system is provided via a function object that evaluates
a set of simultaneous first-order differential equations. The input
to the functions is a pair that contains the independent parameter and
corresponding collection of dependent parameter values.

Compatible equation systems may be implemented by inheriting from the
abstract baseclass, aply::math::DiffEqSystem.

\par Example
\snippet test/test_DiffEqSolve.cpp DoxyExample00

*/

class DiffEqSolve
{

private: // disable

	//! Disable implicit copy
	DiffEqSolve
		(DiffEqSolve const & orig);

	//! Disable implicit assignment
	DiffEqSolve &
	operator=
		(DiffEqSolve const & rhs);

public: // methods
	//! Construct with a given stepsize.
	explicit
	DiffEqSolve
		(double const & stepSize);

	std::pair<double, std::vector<double> >
	solutionFor
		( double const & xValue
		, math::DiffEqSystem const & equations
 		) const;

	// copy constructor -- compiler provided
	// assignment operator -- compiler provided
    // destructor -- compiler provided

	//! Descriptive information about this instance.
	std::string
	infoString
		( std::string const & title=std::string()
		, std::string const & fmt="%20.15g"
		) const;

private:

	int
	rk4
	( double const & stop 
	, math::DiffEqSystem const & functor
	, std::pair<double, std::vector<double> > & result
	) const;

	double theStep;
};

} // [math]
} // [aply]

#endif // aply_math_DiffEqSolve_INCL_
