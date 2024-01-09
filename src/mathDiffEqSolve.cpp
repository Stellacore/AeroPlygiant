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

/*! \file
\brief Definitions for math::DiffEqSolve
*/


#include "mathDiffEqSolve.hpp"

#include <Engabra>

#include <algorithm>
#include <sstream>


namespace aply
{
namespace math
{

DiffEqSolve :: DiffEqSolve
	( double const & stepSize
	)
	: theStep(stepSize)
{
}

std::pair<double, std::vector<double> >
DiffEqSolve :: solutionFor
	( double const & xValue
	, DiffEqSystem const & equations
	) const
{
	std::pair<double, std::vector<double> > const initValues
		(equations.initValues());
	std::pair<double, std::vector<double> > out(initValues);

	int istat = rk4 (xValue, equations, out);
	if (istat != 0)
	{
		std::cerr << "DiffEqSolve::eval: rk4() error status = "
			<< istat << std::endl;
		return initValues;
	}

	return out;
}

int
DiffEqSolve :: rk4
	( double const & stop
	, DiffEqSystem const & functor
	, std::pair<double, std::vector<double> > & result
	) const
{
	std::pair<double, std::vector<double> > const init
		(functor.initValues());
	
	double const start(init.first);
	double const step(theStep);

	std::vector<double> yVec(init.second);

	std::vector<double> F1;
	std::vector<double> F2;
	std::vector<double> F3;
	std::vector<double> F4;

	std::vector<double> tmp(yVec.size());

//
// XXX -- separate stepping control from RK formulae
//

	bool done(false);
	int nstep(0);
	double delta(std::abs(step));
	if (stop < start)
	{ 
		delta = -delta;
	}
	double delo2(delta/2.0);
	double delo6(delta/6.0);
	while (! done)
	{
// XXX needs comments to clarify special case
		double const tparm(start + nstep * delta);
		if (std::abs(stop-tparm) < std::abs(delta))
		{
			delta = stop-tparm;
			delo2 = delta/2.0;
			delo6 = delta/6.0;
			done = true;
		}

//
// XXX -- double check: looks like the abssissa update was missing
//        e.g. (tparm + XXX) for K2,3,4
//

		// RK coefficients

		// K1 = f(xn, yn)
		F1 = functor(std::make_pair(tparm, yVec));

// XXX BTW if efficiency were a concern here (which it is not), the
//     multiply/add combo would be more effectively evaluated with
//     something like an xstl::addMultiple() transform that does both the
//     add and multiply with only 1x the loop overhead. Efficiency
//     would be a bad argument in this case... but clarity of code
//     might be a pretty good one... e.g.
//
// xstl::addMultiple(yVec.begin(), yVec.end(), del02, F1.begin(), out.begin());

		// K2 = f(xn+h/2 , yn+K1*h/2)
		/*
		xstl::multiply(F1.begin(), F1.end(), delo2, tmp.begin());
		xstl::add(yVec.begin(), yVec.end(), tmp.begin(), tmp.begin());
		*/

		std::transform
			( F1.begin(), F1.end()
			, tmp.begin()
		//	, std::bind2nd(std::multiplies<double>(), delo2) );
			, [delo2] (double const & val) { return delo2 * val; }
			);

		std::transform
			( yVec.begin(), yVec.end()
			, tmp.begin()
			, tmp.begin()
			, std::plus<double>() );

		F2 = functor(std::make_pair(tparm + delo2, tmp));

		// K3 = f(xn+h/2 , yn+K2*h/2)
		/*
		xstl::multiply(F2.begin(), F2.end(), delo2, tmp.begin());
		xstl::add(yVec.begin(), yVec.end(), tmp.begin(), tmp.begin());
		*/
		std::transform
			( F2.begin(), F2.end()
			, tmp.begin()
		//	, std::bind2nd(std::multiplies<double>(), delo2) );
			, [delo2] (double const & val) { return delo2 * val; }
			);

		std::transform
			( yVec.begin(), yVec.end()
			, tmp.begin()
			, tmp.begin()
			, std::plus<double>() );

		F3 = functor(std::make_pair(tparm + delo2, tmp));

		// K4 = f(xn+h , yn+K3*h)
		/*
		xstl::multiply(F3.begin(), F3.end(), delta, tmp.begin());
		xstl::add(yVec.begin(), yVec.end(), tmp.begin(), tmp.begin());
		*/
		std::transform
			( F3.begin(), F3.end()
			, tmp.begin()
		//	, std::bind2nd(std::multiplies<double>(), delta) );
			, [delta] (double const & val) { return delta * val; }
			);

		std::transform
			( yVec.begin(), yVec.end()
			, tmp.begin()
			, tmp.begin()
			, std::plus<double>() );

		F4 = functor(std::make_pair(tparm + delta, tmp));

		// Solution update

		// 2.*K2 + 2.*K3
		/*
		xstl::add(F2.begin(), F2.end(), F3.begin(), tmp.begin());
		xstl::multiply(tmp.begin(), tmp.end(), 2.0, tmp.begin());
		*/
		std::transform
			( F2.begin(), F2.end()
			, F3.begin()
			, tmp.begin()
			, std::plus<double>() );

		std::transform
			( tmp.begin(), tmp.end()
			, tmp.begin()
		//	, std::bind2nd(std::multiplies<double>(), 2.) );
			, [] (double const & val) { return 2. * val; }
			);

		// F1 + (2.*K2 + 2.*K3)
		//xstl::add(F1.begin(), F1.end(), tmp.begin(), tmp.begin());
		std::transform
			( F1.begin(), F1.end()
			, tmp.begin()
			, tmp.begin()
			, std::plus<double>() );

		// (F1 + 2.*K2 + 2.*K3) + F4
		//xstl::add(tmp.begin(), tmp.end(), F4.begin(), tmp.begin());
		std::transform
			( tmp.begin(), tmp.end()
			, F4.begin()
			, tmp.begin()
			, std::plus<double>() );

		// (h/6) * (F1 + 2.*K2 + 2.*K3 + F4)
		//xstl::multiply(tmp.begin(), tmp.end(), delo6, tmp.begin());
		std::transform
			( tmp.begin(), tmp.end()
			, tmp.begin()
		//	, std::bind2nd(std::multiplies<double>(), delo6) );
			, [delo6] (double const & val) { return delo6 * val; }
			);

		// yn + (h/6) * (F1 + 2.*K2 + 2.*K3 + F4)
		//xstl::add(yVec.begin(), yVec.end(), tmp.begin(), yVec.begin());
		std::transform
			( yVec.begin(), yVec.end()
			, tmp.begin()
			, yVec.begin()
			, std::plus<double>() );

		++nstep;
	}

	result = make_pair(stop, yVec);

	return 0;
}

//
// infoString
//
std::string
DiffEqSolve :: infoString
	( std::string const & title
	, std::string const & fmt
	) const
{
	std::ostringstream oss;

	if (! title.empty())
	{
		oss << title << std::endl;
	}

	oss << "Step size: " << engabra::g3::io::fixed(theStep);
	return oss.str();
}

} // [math]
} // [aply]
