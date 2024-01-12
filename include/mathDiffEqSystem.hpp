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

#ifndef aply_math_DiffEqSystem_INCL_
#define aply_math_DiffEqSystem_INCL_


/*! \file
\brief Declarations for math::DiffEqSystem
*/


#include <string>
#include <utility>
#include <vector>


namespace aply
{
namespace math
{

	/* TODO - maybe useful idea, but needs explained better and demonstrated.
	To include the independent parameter in equation system, introduce a 
	function of this variable and insert it into the simultaneous system.
	x' = 1.

	If you need to include x itself in the equations, you can make it into a
	a dependent variable x' = 1. For the initial value, you should repeat the 
	value of the independant parameter as the value of x.
	*/

/*! \brief Abstract class for system of first order differential equations.

\b Summary:

In general, the functor, operator(), must transform a system of function
values (y1, y2, y3, ...) into a system of derivative values (y1', y2',
y3', ...).  The initValues() function must be provide values for the
functions (y1, y2, y3, ...) that are associated with a specific initial
value for the (assumed) independent parameters, x.

Once a desired ODE System is expressed within a derived class, then the
system may be solved using aply::math::DiffEqSolve.

\b Example:

Ref: Concrete class aply::examp::diffeq::UniformAccel
in example/diffeqSystem.hpp which implements a uniform acceleration
as a scalar valued equation system (y''=const).

The ODE system is implemented as:
\snippet example/diffeqSystem.hpp DoxyExampleDiffEqSystemOp()

The initial value function is implemented as:
\snippet example/diffeqSystem.hpp DoxyExampleDiffEqSystemInitVal

\b Detail:

Let x be an independent parameter (e.g. often time or distance). Also
let y=y(x) be the function of interest (i.e. the solution function). Let
y'(x), y''(x), ... be increasingly higher order derivatives (with respect
to x) of the function y(x).

A typical ordinary differential equation (ODE) problem involves a single
(vector or scalar valued) equation that expresses the relationship
between derivatives of various orders. E.g. an equation of the form
f(x,y,y',...,y[n]) represents an n-th order equation.

For this algorithm the function must be explicit in the n-th order term.
I.e. f() can be rearranged to be expressed as y[n]=g(x,y,y'...,y[n-1]).

An equation like this can be expressed as a system of simultaneous
first order ODE equations:
\code
y[n] = g(x,y,y',...,y[n-1])
y[n-1] = (y[n-2])'
y[n-2] = (y[n-1])'
...
y[2] = y'' = (y')' = y[1]'
y[1] = y'  = (y)'  = desired answer function
\endcode


If you have an ODE y' = x^2, add x' = 1 to this to make a system of two
equations. From the functor you would return (1.0, in[0] * in[0]) if you put
x first.

\par Scalar first-order equation
A single, first order differential equation is represented as
\code
y' = f(x,y)
\endcode
where f(x) is the known derivative function that is to be integrated.

\par Vector first-order equation
A vector equation simply repeats the scalar equation structure for
each of the vector component. E.g. consider finding the position in 3D of
a particle as a function of time, x(t), when the velocity function, v(t)
is known.
\code
v1(t) = x1' = f1(t,v1,v2,v3)
v2(t) = x2' = f2(t,v1,v2,v3)
v3(t) = x3' = f3(t,v1,v2,v3)
\endcode
Note that the function for the velocity of each component (e.g. vi(t)) can
depend on the velocity of the other components. So that, in general, these
equations are coupled.

\par Scalar second-order equation
A single, second-order differential equation may be decomposed into two
coupled first-order equations. e.g. For example, consider the second-order
equation
\code
y'' = g(t,y,y')
where y=y(t) is the desired solution
\endcode
Let
\code
y0  = f0(t) is the desired solution
y0' = y1
y1' = g(t,y0,y1)
\endcode
This system of equations may be expressed as
\code
y0  = f0(t, y0, y1) = solution
y0' = f1(t, y0, y1) = y1
y1' = f2(t, y0, y1) = g(t,y0,y1)
\endcode

\par Vector higher-order equations
The process illustrated above for reducing a second-order equation to
coupled first order equations can be continued to any order. Furthermore,
the process can be combined with decomposition of vector equations into
coupled component equations as illustrated above. Together, these
processes can be applied repeatedly to reduce an arbitrary high-order
vector differential equation to a set of coupled simultaneous first-order
scalar equations.

\par Example of Vector-valued 2nd order ODE system:

Functor implementation:
\snippet demo/demoIntegrate.cpp DoxyExample00

Initial values implementation:
\snippet demo/demoIntegrate.cpp DoxyExample01


*/
/* - no test

\par Example
\dontinclude test/test_DiffEqSystem.cpp
\skip ExampleStart
\until ExampleEnd
 */
class DiffEqSystem
{
public:

	//! typical dtor.
	virtual
	~DiffEqSystem
		()
	{ }

	/*! \brief Values of all derivative functions evaluated at xyValues.
	 *
	 * The return vector contains the values of the first order
	 * derivatives of functions: f0, f1, f2, .... I.e.
	 * \arg return std::vector<double>{ f0', f1', f2', ... };
	 *
	 * \par Example
	 * Consider the pure 3rd order differential equation
	 * \arg y''' = f3(x, y, y', y'')
	 *
	 * The input arguments in #xyValues are
	 * \arg x   = xyValues.first
	 * \arg y   = f0 = xyValues.second[0]
	 * \arg y'  = f1 = xyValues.second[1]
	 * \arg y'' = f2 = xyValues.second[2]
	 *
	 * The return values are the first order derivatives, e.g.
	 * \arg f1 = f0'
	 * \arg f2 = f1'
	 * \arg f3 = f3(x, f0, f1, f2)
	 *
	 * or in code form, should then look like
	 * \code
	 * return std::vector<double>
	 * 	{ f1
	 * 	, f2
	 * 	, f3 = f3(x, f0, f1, f2)
	 *	};
	 * \endcode
	 */
	virtual
	std::vector<double>
	operator()
		( std::pair<double, std::vector<double> > const & xyValues
	 		//!< values represent [(x), (y0, y1, ...)]
		) const = 0;

	/*! \brief Initial conditions for the system of equations.
	 *
	 * The return values represent: [(x), (y0, y1, ...)]
	 *
	 * For the example in operator() description, code would be  like
	 * \code
	 * return std::vector<double>
	 * 	{ t0
	 * 	, f0(t0)
	 * 	, f1(t0)
	 *	};
	 * \endcode
	 */
	virtual
	std::pair< double, std::vector<double> >
	initValues
		() const = 0;
};

} // [math]
} // [aply]

#endif // aply_math_DiffEqSystem_INCL_
