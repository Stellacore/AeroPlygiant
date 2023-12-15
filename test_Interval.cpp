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
 * \brief Unit test for class geom::Interval
 *
 */


#include "tst.hpp"

#include "geomInterval.hpp"

#include <sstream>


namespace
{
	//! Check .... TODO
	void
	test0
		( std::ostringstream & oss
		)
	{
		// [DoxyExample00]
		// example (half-open) interval
		aply::geom::Interval gap23(2., 3.);

		// example values
		double const gotValueAt1{ gap23.fracAtValue( 1.) };
		double const gotValueAt2{ gap23.fracAtValue( 2.) };
		double const gotValueAt3{ gap23.fracAtValue( 3.) };
		double const gotValueAt4{ gap23.fracAtValue( 4.) };

		// associated fractions
		double const expValueAt1{ -1. }; // extrapolation before interval
		double const expValueAt2{  0. }; // start point IN-cluded
		double const expValueAt3{  1. }; // end point EX-cluded
		double const expValueAt4{  2. }; // extrapolation after interval
		// [DoxyExample00]

		using tst::checkGotExp;

		checkGotExp(oss, gotValueAt1, expValueAt1, "fracAtValue");
		checkGotExp(oss, gotValueAt2, expValueAt2, "fracAtValue");
		checkGotExp(oss, gotValueAt3, expValueAt3, "fracAtValue");
		checkGotExp(oss, gotValueAt4, expValueAt4, "fracAtValue");

		double const rtnValueAt1{ gap23.valueAtFrac(expValueAt1) };
		double const rtnValueAt2{ gap23.valueAtFrac(expValueAt2) };
		double const rtnValueAt3{ gap23.valueAtFrac(expValueAt3) };
		double const rtnValueAt4{ gap23.valueAtFrac(expValueAt4) };

		checkGotExp(oss, rtnValueAt1, 1., "valueAtFrac");
		checkGotExp(oss, rtnValueAt2, 2., "valueAtFrac");
		checkGotExp(oss, rtnValueAt3, 3., "valueAtFrac");
		checkGotExp(oss, rtnValueAt4, 4., "valueAtFrac");

	}
}


/*! \brief Unit test for geom::Interval
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);

	return tst::finish(oss);
}

