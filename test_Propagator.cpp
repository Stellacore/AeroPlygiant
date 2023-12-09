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
 * \brief Unit test for class Propagator
 *
 */


#include "env.hpp"
#include "ray.hpp"

#include "tst.hpp"


namespace
{
	//! True if (minIncluded <= value < maxExcluded)
	inline
	bool
	inInterval
		( double const & minIncluded
		, double const & value
		, double const & maxExcluded
		)
	{
		return ( (! (value < minIncluded)) && (value < maxExcluded) );
	}

	using namespace engabra::g3;

	//! Simple test volume
	struct UnitBox : public env::ActiveVolume
	{
		virtual
		bool
		contains
			( Vector const & rVec
			) const
		{
			return
				(  inInterval(0., rVec[0], 1.)
				&& inInterval(0., rVec[1], 1.)
				&& inInterval(0., rVec[2], 1.)
				);
		}

	}; // UnitBox

	// Construct media (providing values outside active volume)
	struct AirCube : public env::IndexVolume
	{
		AirCube
			()
			: IndexVolume(UnitBox{})
		{}

		//! Active volume should restrict use of indices
		virtual
		double
		nuValue
			( Vector const & rVec
			) const
		{
			constexpr double nuAirSTP{ 1.000273 };
			return nuAirSTP;
		}

	}; // AirCube


	//! Check ray propagation through a small uniform volume
	void
	test0
		( std::ostream & oss
		)
	{

		// [DoxyExample00]
		// [DoxyExample00]

		// construct media environment
		AirCube const opticalMedia; // cube of 'air' (nu=1.000271)

		// configure propagator
		constexpr double propStepSize{ 1./16. };
		ray::Propagator const prop{ &opticalMedia, propStepSize };

		// configure the ray(s) for propagation
		// start ray heading in +x direction, starting at x=-.5
		ray::Start const start{ ray::Start::from(e1, -.5*e1 ) };

		// trace the ray(s)
		constexpr double saveDeltaDistance{ 1./8. };
		ray::Path aPath(start, 2.*e1, saveDeltaDistance);
		prop.traceNodes(start, &aPath);

oss << "Failure\n";
	}
}


/*! \brief Unit test for Propagator
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);

	return tst::finish(oss);
}

