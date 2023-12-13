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

#ifndef aply_geom_Cylinder_INCL_
#define aply_geom_Cylinder_INCL_

/*! \file
 *
 * \brief Geometric utilities.
 *
 */


#include <Engabra>


namespace aply
{
namespace geom
{
	using namespace engabra::g3;

	//! \brief A geometric cylinder shape of finite length
	struct Cylinder
	{
		Vector const theAxisBeg;
		Vector const theAxisDir;
		double const theLength;
		double const theRadius;
		Interval const theGapLength;
		Interval const theGapRad;

		//! Value construction
		inline
		explicit
		Cylinder
			( Vector const & axisBeg
			, Vector const & axisDir
			, double const & length
			, double const & radius
			)
			: theAxisBeg{ axisBeg }
			, theAxisDir{ direction(axisDir) }
			, theLength{ length }
			, theRadius{ radius }
			, theGapLength{ 0., theLength }
			, theGapRad{ 0., theRadius }
		{ }

		//! Distance orthogonal from body axis to loc.
		inline
		double
		distanceFromAxis
			( Vector const & someLoc
			) const
		{
			double const dist
				{ ((someLoc - theAxisBeg)*theAxisDir).theSca[0] };
			return dist;
		}

		//! Distance orthogonal from body axis to loc.
		inline
		double
		fractionFromAxis
			( Vector const & someLoc
			) const
		{
			return theGapRad.fracAtValue(distanceFromAxis(someLoc));
		}

		//! Distance parallel along body axis to loc.
		inline
		double
		distanceAlongAxis
			( Vector const & someLoc
			) const
		{
			double const dist
				{ (someLoc*theAxisDir).theSca[0] };
			return dist;
		}

		inline
		double
		fractionAlongAxis
			( Vector const & someLoc
			) const
		{
			return theGapLength.fracAtValue(distanceAlongAxis(someLoc));
		}

	}; // Cylinder

} // geom

} // [aply]

#endif // aply_geom_Cylinder_INCL_

