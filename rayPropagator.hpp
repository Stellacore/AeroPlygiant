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

#ifndef aply_ray_Propagator_INCL_
#define aply_ray_Propagator_INCL_

/*! \file
 *
 * \brief Ray propagation simulation functions.
 *
 */


#include "rayDirChange.hpp"
#include "rayNode.hpp"

#include "env.hpp"

#include <Engabra>

#include <cmath>
#include <limits>
#include <utility>


/*! \brief Functions and classes for simulation of ray propagation
 */
namespace ray
{
	using namespace engabra::g3;

	//! Update tangent direction across single (idealized) interface boundary.
	inline
	std::pair<Vector, DirChange>
	nextTangentDir
		( Vector const & tDirPrev //!< Must be unit length
		, double const & nuPrev //!< Incoming IoR
		, Vector const & gCurr //! Must be non-zero (to be invertable
		, double const & nuNext //!< Exiting IoR
		)
	{
		// default case is an unaltered ray
		std::pair<Vector, DirChange> tanDirChange{ tDirPrev, Unaltered };
		Vector & tDirNext = tanDirChange.first;
		DirChange & tChange = tanDirChange.second;
		//
		// check for stop condition
		if (! engabra::g3::isValid(nuPrev))
		{
			tChange = Stopped;
		}
		else
		{
			// compute refraction bivector
			// note that magnitude is order of |gCurr|
			BiVector const currB{ (nuPrev/nuNext) * (tDirPrev*gCurr).theBiv };
			//
			// note that sq(bivector) = -magSq(bivector)
			double const gCurrSq{ magSq(gCurr) };
			double const radicand{ gCurrSq - magSq(currB) };
			//
			// use current conditions to select computation option
			//
			Vector const gCurrInv{ (1./gCurrSq) * gCurr };
			if (radicand < 0.) // total internal reflection
			{
				// reflect tangent from interface plane (dual to gCurr)
				tDirNext = -(gCurr * tDirPrev * gCurrInv).theVec;
				tChange = Reflected;
			}
			else
			{
				double const rootXi{ std::sqrt(radicand) };
				double const tDotG{ (tDirPrev * gCurr).theSca[0] };
				if (tDotG < 0.) // propagating into less dense media
				{
					Spinor const spin{ -rootXi, currB };
					tDirNext = (spin * gCurrInv).theVec;
					tChange = Diverged;
				}
				else
				if (0. < tDotG) // propagating into more dense media
				{
					Spinor const spin{  rootXi, currB };
					tDirNext = (spin * gCurrInv).theVec;
					tChange = Converged;
				}
				// (nuNext == nuPrev) // same as default (gCurr == 0)
			}
		}
		//
		return tanDirChange;
	}

	//! Ray propagation functions
	struct Propagator
	{
		env::IndexVolume const * const thePtMedia{ nullptr };
		double const theStepDist{ null<double>() };

	private:

		struct Step // Propagator::
		{
			double theNextNu;
			Vector theNextTan;
			DirChange theChange;

		}; // Step

		//! Estimate next tangent based on local object refraction
		inline
		Step
		nextStep // Propagator::
			( Vector const & tPrev //!< Must be unit length
			, double const & nuPrev
			, Vector const & rCurr
			) const
		{
			double nuNext{ null<double>() };
			Vector tNext{ tPrev }; // iteratively evolved from here
			DirChange change{ Null };
			//
			// Check if there's anything to compute (vs unaltered propagation)
			Vector const gCurr{ thePtMedia->nuGradient(rCurr, theStepDist) };
			double const gMag{ magnitude(gCurr) };
			static double const gTol // enough to unitize and invert gCurr
				{ std::numeric_limits<double>::min() };
			if (! (gTol < gMag)) // unaltered
			{
				// location at which to evaluate nuNext (iteratively refined)
				Vector const qNext{ rCurr + .5*theStepDist*tNext };

				// update estimated forward next refraction index value
				// (which should be same as previous)
				nuNext = thePtMedia->qualifiedNuValue(qNext);

				// update tangent direction
				tNext = tPrev; // default initialized

				// update boundary action
				change = Unaltered; // default initialized
			}
			else
			// if (gTol < gMag) // ray path (tangent direction) changes
			{
				bool isReflection{ false }; // used to exit early 
				//
				// iterate on determination of exit media IoR
				// (ref Refraction.lyx doc)
				//
				// location at which to evaluate nuNext (iteratively refined)
				Vector qNext{ rCurr + .5*theStepDist*tNext };
				//
				// tolerance until epsilon < difSq (sqrt(eps)<|dif|)
				static double const tolDifSq
					{ std::numeric_limits<double>::epsilon() };
				double difSq{ 2.*tolDifSq }; // large value to force first loop
				constexpr std::size_t maxLoop{ 10u }; // avoid infinite loop
				std::size_t numLoop{ 0u };
				bool doLoop{ true };
				while ((tolDifSq < difSq) && (numLoop++ < maxLoop) && doLoop)
				{
					// update IoR evaluation location
					if (! isReflection)
					{
						// update refraction index to midpoint of predicted next
						// interval (along evolving next tangent direction).
						qNext = rCurr + .5*theStepDist*tNext;
					}
					else
					// if (isReflection) // perfect reflection
					{
						Vector const gDir{ direction(gCurr) };
						qNext = rCurr + .5*theStepDist*gDir;
						// No need to iterate further for perfect reflection
						doLoop = false; // exit after updating next ray step
					}

					// update estimated forward next refraction index value
					nuNext = thePtMedia->qualifiedNuValue(qNext);

					// update tangent direction
					std::pair<Vector, DirChange> const tDirChange
						{ nextTangentDir(tPrev, nuPrev, gCurr, nuNext) };

					// check for stop condition
					if (Stopped == tDirChange.second)
					{
						break;
					}

					// note reflection condition for next loop iteration
					isReflection = (Reflected == tDirChange.second);

					// evaluate convergence of tangent direction
					Vector const & tResult = tDirChange.first;
					difSq = magSq(tResult - tNext);

					// update tangent direction
					tNext = tResult;
					change = tDirChange.second;

				} // while loop on refraction index estimation

			} // significant gCurr magnitude

			// check for invalid media index volume (e.g. exit region)
			if (! engabra::g3::isValid(nuNext))
			{
				change = Stopped;
			}

			return Step{ nuNext, tNext, change };
		}

		//! Predicted next location stepsize units along tangent from rVec
		inline
		Vector
		nextLocation // Propagator::
			( Vector const & rVec
			, Vector const & tVec
			) const
		{
			return { rVec + theStepDist * tVec };
		}

	public:

		//! True if this instance is valid
		inline
		bool
		isValid // Propagator::
			() const
		{
			return engabra::g3::isValid(theStepDist);
		}

		/*! Perform forward integration step by step.
		 *
		 * Essentially is Euler's method for integration of the ray path
		 * (with all attendent pitfalls).
		 *
		 */
		template <typename Consumer>
		inline
		void
		tracePath // Propagator::
			( Consumer * const & ptConsumer
			) const
		{
			if (isValid() && ptConsumer)
			{
				Vector const & tBeg = ptConsumer->theStart.theTanDir;
				Vector const & rBeg = ptConsumer->theStart.thePntLoc;

				// start with initial conditions
				Vector tPrev{ tBeg };
				Vector rCurr{ rBeg };

				// incident media IoR
				Vector const rPrev{ (rBeg - .5*theStepDist*tBeg) };
				double nuPrev{ thePtMedia->qualifiedNuValue(rPrev) };

				// propagate until path approximate reaches requested length
				// or encounteres a NaN value for index of refraction
				while (ptConsumer->size() < ptConsumer->capacity())
				{
					// determine propagation change at this step
					Step const stepNext{ nextStep(tPrev, nuPrev, rCurr) };
					Vector const & tNext = stepNext.theNextTan;
					double const & nuNext = stepNext.theNextNu;
					DirChange const & change = stepNext.theChange;

					// check for ray termination condition
					if (Stopped == change)
					{
						break;
					}

					// propagate ray to next node location
					Vector const rNext{ nextLocation(rCurr, tNext) };

					// give consumer opportunity to record node data
					Node const nextNode
						{ tPrev, nuPrev, rCurr, nuNext, tNext, change };
					ptConsumer->emplace_back(nextNode);

					// update state for next node
					tPrev = tNext;
					rCurr = rNext;
					nuPrev = nuNext;
				}
			}
		}

	}; // Propagator

} // [ray]



#endif // aply_ray_Propagator_INCL_

