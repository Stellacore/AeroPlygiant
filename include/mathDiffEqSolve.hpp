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

/*! \brief This class uses numerical methods to solve ordinary
differential equations.

The function is provided via a function object that evaluates a set of
simultaneous first-order differential equations. The input to the functions
is a pair that contains the independent parameter and corresponding collection
of dependent parameter values.

\par Example
\dontinclude testmath/uDiffEqSolve.cpp
\skip ExampleStart
\until ExampleEnd

XXX  -- update comments from here

Here note that you must specify the step size used by the iterative method.
Additionally, you must specify initial values to fix the solution.

You must write a functor to specify the differential equation. The
implied prototype of the functor is std::vector<double>
functor(std::vector<double> const & in). For a first-order ODE
that does not use x (the independent parameter), the input vector is of size 1
and contains the value of y. The functor must return a vector of size 1
containing the value of y'. So for the differential equation y' = y, you
would return a vector ( in[0] ).

You also must use the init function to specify the value of the function at
a point to fix the solution. Here, you could use the point 0.0 and require
the value to be 1.0, which will fix the function as exp(x).

For higher-order ODEs (again assuming x is not used), you must rewrite the
ODE as a system of first-order ODEs. So, given y'' = -y, rewrite this as
y' = z and z' = -y. Any system can be straightforwardly converted in this
fashion.

Now the functor will be passed a vector (y, z) of the two function values.
You must return the vector (y', z'). For this example, you would return
(-in[1], in[0] ). For the initial values, you specify the values (y, z) at a
given point. So you could use the point 0.0 and pass in (0.0, 1.0). This
specifies that y(0) = 0 and that z(0) = 1. Since z = y', this means that
y'(0) = 1. This fixes the function as sin(x).

To vary the phase and amplitude, you can return (0.0, amplitude) using a
point "phase" where you want the curve to be 0. This will specify the function
as amplitude * sin(x - phase).

In general, the functor must transform (y1, y2, y3, ...) into
(y1', y2', y3', ...) and the init function must be initialized with values
(y1, y2, y3, ...) at the point x.

If you need to include x itself in the equations, you can make it into a
a dependent variable x' = 1. For the initial value, you should repeat the 
value of the independant parameter as the value of x.

If you have an ODE y' = x^2, add x' = 1 to this to make a system of two
equations. From the functor you would return (1.0, in[0] * in[0]) if you put
x first.
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
