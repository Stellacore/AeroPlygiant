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
		explicit
		UnitBox
			()
			: ActiveVolume("UnitBox")
		{ }

		//! Create box with unit dimensions
		inline
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

	static std::shared_ptr<env::ActiveVolume>
		const sPtUnitBox{ std::make_shared<UnitBox>() };

	// Construct media (providing values outside active volume)
	struct AirCube : public env::IndexVolume
	{
		AirCube
			()
			: IndexVolume(sPtUnitBox)
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
	testBox
		( std::ostream & oss
		)
	{
		UnitBox const box;
		std::vector<Vector> const inLocs
			{ Vector{ 0., 0., 0. } // start corner is in
			, Vector{ .5, .5, .5 } // interior point is in
			};
		std::vector<Vector> const outLocs
			{ Vector{ 1., 1., 1. } // end corner is out
			, Vector{-.5,-.5,-.5 } // exterior point is out
			, Vector{1.5,1.5,1.5 } // exterior point is out
			};
		for (Vector const & inLoc : inLocs)
		{
			if (! box.contains(inLoc))
			{
				oss << "Failure of inloc test\n";
				oss << "inLoc: " << inLoc << '\n';
			}
		}
		for (Vector const & outLoc : outLocs)
		{
			if (  box.contains(outLoc))
			{
				oss << "Failure of outloc test\n";
				oss << "outLoc: " << outLoc << '\n';
			}
		}
	}

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
		constexpr double propStepSize{ 1./8. };
		ray::Propagator const prop{ &opticalMedia, propStepSize };

		// configure the ray(s) for propagation
		// start ray heading in +x direction, starting at x=0.
		ray::Start const start{ ray::Start::from(e1, zero<Vector>() ) };

		// trace the ray(s)
		constexpr double saveDeltaDistance{ 1./8. };
Vector const stopLoc{ 2. * e1 };
//Vector const stopLoc{ zero<Vector>() };
//Vector const stopLoc{ null<Vector>() };
		ray::Path aPath(start, stopLoc, saveDeltaDistance);
		prop.tracePath(&aPath);

		std::cout << aPath.infoString("aPath") << std::endl;

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

	testBox(oss);
	test0(oss);

	return tst::finish(oss);
}

