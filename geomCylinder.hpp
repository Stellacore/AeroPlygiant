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


#include "geomInterval.hpp"

#include <Engabra>


namespace aply
{
namespace geom
{
	using namespace engabra::g3;

	//! \brief A geometric cylinder shape of finite length
	struct Cylinder
	{
		//! Start point of axis.
		Vector const theAxisBeg;
		//! (unitary) direction of axis leaving #theAxisBeg.
		Vector const theAxisDir;
		//! Length of cylinder along the axis (distance between end caps).
		double const theLength;
		//! Radius of cylinder (axis to outer curved edge).
		double const theRadius;
		//! Length Interval from begin cap to end cap.
		Interval const theLengthInterval;
		//! Radial Interval from axis to outer curved edge.
		Interval const theRadialInterval;

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
			, theLengthInterval{ 0., theLength }
			, theRadialInterval{ 0., theRadius }
		{ }

		//! Distance orthogonal from body axis to loc.
		inline
		double
		distanceFromAxis
			( Vector const & someLoc
			) const
		{
			Vector const relLoc{ someLoc - theAxisBeg };
			// dot product of relLoc and axis direction
			BiVector const rejection{ (relLoc*theAxisDir).theBiv };
			return magnitude(rejection);
		}

		//! Distance orthogonal from body axis to loc.
		inline
		double
		fractionFromAxis
			( Vector const & someLoc
			) const
		{
			return theRadialInterval.fracAtValue(distanceFromAxis(someLoc));
		}

		//! Distance parallel along body axis to loc.
		inline
		double
		distanceAlongAxis
			( Vector const & someLoc
			) const
		{
			Vector const relLoc{ someLoc - theAxisBeg };
			// dot product of relLoc and axis direction
			double const projection{ (relLoc*theAxisDir).theSca[0] };
			return projection;
		}

		inline
		double
		fractionAlongAxis
			( Vector const & someLoc
			) const
		{
			return theLengthInterval.fracAtValue(distanceAlongAxis(someLoc));
		}

	}; // Cylinder

} // geom

} // [aply]

#endif // aply_geom_Cylinder_INCL_

