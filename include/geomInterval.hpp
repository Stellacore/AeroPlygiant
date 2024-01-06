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

#ifndef aply_geom_Interval_INCL_
#define aply_geom_Interval_INCL_

/*! \file
 *
 * \brief Geometric utilities.
 *
 */


#include <utility>


namespace aply
{
namespace geom
{

	/*! \brief Use to values to define a distance scale (origin and unit value)
	 *
	 * Perhaps best explained by example:
	 * \snippet test_Interval.cpp DoxyExample00
	 *
	 * And the inverse, e.g.
	 * \arg valueAtFrac(.75, pair<>(2., 3.)) == 2.75
	 *
	 * \note The "include/exclude" conditions are not relevant to
	 *       currently implemented methods.
	 *
	 */
	class Interval
	{
		//! Define the half open interval [min,max)
		std::pair<double, double> const theMinMax{};
		//! Distance between end points.
		double const theSpan{};
		//! Inverse of the span (i.e. scale = 1./span).
		double const theScale{};

	public:

		//! Value construction of half open interval [begValue, endValue)
		explicit
		Interval
			( double const & begValue //!< IN-cluded value starting interval
			, double const & endValue //!< EX-cluded value ending interval
			)
			: theMinMax{ begValue, endValue }
			, theSpan{ theMinMax.second - theMinMax.first }
			, theScale{ 1. / theSpan }
		{ }

		//! The origin of the interval
		inline
		double
		min
			() const
		{
			return theMinMax.first;
		}

		//! The origin of the interval
		inline
		double
		max
			() const
		{
			return theMinMax.second;
		}

		/*! \brief Inter(extra)polated fraction of way into interval.
		 */
		inline
		double
		fracAtValue
			( double const & value
			) const
		{
			return theScale * (value - min());
		}

		//! \brief Value associated with fraction between end points.
		inline
		double
		valueAtFrac
			( double const & frac
			) const
		{
			return (frac * theSpan + min());
		}

	}; // Interval

} // geom

} // [aply]

#endif // aply_geom_Interval_INCL_

