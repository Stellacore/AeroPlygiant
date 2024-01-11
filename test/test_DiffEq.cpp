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
 * \brief Test validity of differential equation formulation.
 *
 */


#include "envActiveVolume.hpp"
#include "envIndexVolume.hpp"
#include "math.hpp"
#include "rayPath.hpp"
#include "rayPropagator.hpp"
#include "rayStart.hpp"

#include "tst.hpp"

#include <Engabra>

#include <cmath>
#include <sstream>
#include <vector>


namespace
{
	//! Vector-Vector exterior product.
	engabra::g3::BiVector
	wedge
		( engabra::g3::Vector const & vecA
		, engabra::g3::Vector const & vecB
		)
	{
		return (vecA * vecB).theBiv;
	}

} // [annon]

namespace tst
{
	using namespace engabra::g3;

	//! \brief A 3D Gaussian concentration
	struct Blob
	{
		Vector const theCenter{ null<Vector>() };
		double const theArgScl{ null<double> () };
		double const theAmp{ null<double> () };

		//! Value construction
		inline
		explicit
		Blob
			( Vector const & center
			, double const & sigma = 1.
			, double const & maxValue = 1.
			)
			: theCenter{ center }
			, theArgScl{ -1. / (2. * sq(sigma)) }
			// , theAmp{ 1. / cube(sigma * std::sqrt(2. * pi)) } // unit volume
			, theAmp{ maxValue } // unit value at center
		{ }

		//! \brief Value of the 3D distribution at rLoc
		inline
		double
		operator()
			( Vector const & rLoc
			) const
		{
			return theAmp * std::exp(theArgScl * magSq(rLoc - theCenter));
		}

	}; // Blob

	/*! \brief Several 3D Guassian concentrations of IoR value.
	 */
	struct MediaBlobs : public aply::env::IndexVolume
	{
		//! Attach env::ActiveVolume to define boundaries of the IoR field.
		explicit
		MediaBlobs
			()
			: IndexVolume
				( std::make_shared<aply::env::ActiveBox>
					( engabra::g3::Vector{  0.,  -5.,  -5. }
					, engabra::g3::Vector{  5.,   5.,   5. }
					)
				)
		{ }

		/*! \brief An IoR field with several 3D Gaussian concentrations
		 */
		inline
		double
		nuValue
			( engabra::g3::Vector const & rVec
			) const
		{
			double rfy{ .000 }; // refractivity (sum over blobs)
			// a bunch of blobs with 'glass-like' *refractivity* near centers
			std::vector<tst::Blob> const blobs
				{ tst::Blob(Vector{ 0.,  1.,  .00 }, 1.,  .500)
				, tst::Blob(Vector{ 2.,  0.,  .75 }, 1.,  .500)
				, tst::Blob(Vector{ 4., -0., -.50 }, 1.,  .500)
				};
			for (tst::Blob const & blob : blobs)
			{
				rfy += blob(rVec);
			}
			return (1.000 + rfy);
		}

		//! Functor access to nuValue
		inline
		double
		operator()
			( engabra::g3::Vector const & rVec
			) const
		{
			return nuValue(rVec);
		}

		//! Hessian matrix numerical approximation at rVec
		inline
		aply::math::Matrix
		nuHessian
			( Vector const & rVec
			, double const relStepDist
			) const
		{
			return aply::math::hessianOf(*this, rVec, relStepDist);
		}

	}; // MediaBlobs

} // [tst]

namespace
{
	// Configure propagation step size and specify path save interval
	constexpr double sPropStepDist{ 1./1000. };
	constexpr double sSaveStepDist{ 1./10. };

	//! Numerically generated path
	std::shared_ptr<aply::ray::Path>
	numericalPath
		( std::shared_ptr<tst::MediaBlobs> const & ptrMedia
		, std::size_t const & numNodes = 64u
		)
	{
		std::shared_ptr<aply::ray::Path> ptrPath{ nullptr };

		using namespace aply;
		using namespace engabra::g3;

		// Initial conditions
		ray::Start const start{ ray::Start::from(e1, zero<Vector>()) };

		// Allocate path
		ptrPath = std::make_shared<ray::Path>(start, sSaveStepDist);
		ptrPath->reserve(numNodes);

		// Trace path through a 3D refraction environment
		ray::Propagator const prop{ ptrMedia.get(), sPropStepDist };
		prop.tracePath(ptrPath.get());

		return ptrPath;
	}

	using namespace engabra::g3;

	//! Numerical estimate of derivative values at a node
	struct NodeDiff
	{
		Vector const theLoc{ null<Vector>() };

		double const theNuVal{ null<double>() };
		double const theNuDot{ null<double>() };

		Vector const theGradVal{ null<Vector>() };
		Vector const theGradDot{ null<Vector>() };

		Vector const theNormVal{ null<Vector>() };
		Vector const theNormDot{ null<Vector>() };

		Vector const theTanVal{ null<Vector>() };
		Vector const theTanDot{ null<Vector>() };

		//! Estimate differential values using numerical differencing
		static
		NodeDiff
		from
			( aply::ray::Node const & node
			, std::shared_ptr<tst::MediaBlobs> const & ptrMedia
			, double const & stepDist // path length difference
			)
		{
			// current node information forming basis for central differences
			Vector const & tanPrev = node.thePrevTan;
			Vector const & locCurr = node.theCurrLoc;
			Vector const & tanNext = node.theNextTan;

			// compute evaluation locations for performing differencing
			Vector const locPrev{ locCurr - .5*stepDist*tanPrev };
			Vector const locNext{ locCurr + .5*stepDist*tanNext };

			// estimate derivative (wrt path length) of IoR field scalar
			double const nuPrev{ ptrMedia->nuValue(locPrev) };
			double const nuNext{ ptrMedia->nuValue(locNext) };

			// estimate derivative (wrt path length) of IoR gradient vector
			Vector const gradPrev
				{ ptrMedia->nuGradient(locPrev, .5*stepDist) };
			Vector const gradNext
				{ ptrMedia->nuGradient(locNext, .5*stepDist) };

			// estimate the change in *direction* of the gradient
			Vector const normPrev{ direction(gradPrev) };
			Vector const normNext{ direction(gradNext) };

			double const scale{ 1. / stepDist };

			return
				{ locCurr
				, .5*(  nuNext +   nuPrev) , scale*(  nuNext -   nuPrev)
				, .5*(gradNext + gradPrev) , scale*(gradNext - gradPrev)
				, .5*(normNext + normPrev) , scale*(normNext - normPrev)
				, .5*( tanNext +  tanPrev) , scale*( tanNext -  tanPrev)
				};
		}

	}; // NodeDiff


	//! Evaluate differential equation forms
	struct DiffEq
	{
		std::shared_ptr<tst::MediaBlobs> const thePtrMedia{ nullptr };

		engabra::g3::BiVector
		operator()
			( aply::ray::Node const & node
			) const
		{
			BiVector value{ null<BiVector>() };

			double const dL{ sPropStepDist };

			NodeDiff const nodeDiff{ NodeDiff::from(node, thePtrMedia, dL) };

			using namespace engabra::g3;

			Vector const & gVal = nodeDiff.theGradVal;
			Vector const & gDot = nodeDiff.theGradDot;
			Vector const & uVal = nodeDiff.theNormVal;
			Vector const & uDot = nodeDiff.theNormDot;
			Vector const & tDot = nodeDiff.theTanDot;
			Vector const & tVal = nodeDiff.theTanVal;
			double const & nuVal = nodeDiff.theNuVal;
			double const & nuDot = nodeDiff.theNuDot;

			double const & nuPrev = node.thePrevNu;
			Vector const & tanPrev = node.thePrevTan;
			Vector const & currLoc = node.theCurrLoc;
			double const & nuNext = node.theNextNu;
			Vector const & tanNext = node.theNextTan;

			{
			BiVector const biv1{ nuDot * wedge(tVal, uVal) };
			BiVector const biv2{ nuVal * wedge(tDot, uVal) };
			value = biv1 + biv2;
std::cout << "biv1: " << io::fixed(biv1) << '\n';
std::cout << "biv2: " << io::fixed(biv2) << '\n';
return value; // TODO complete test
			}
			/*
			*/

			// okay
			// Finite difference definition used in path propagation
			BiVector const biv1{ nuPrev * (tanPrev * gVal).theBiv };
			BiVector const biv2{ nuNext * (tanNext * gVal).theBiv };
			value = biv1 - biv2;
std::cout << "biv1: " << io::fixed(biv1) << '\n';
std::cout << "biv2: " << io::fixed(biv2) << '\n';

return value; // TODO complete test

			// Evaluate with finite difference relationships
			double const nuDelta{ nuNext - nuPrev };
			Vector const tanDelta{ tanNext - tanPrev };
			Vector const normDelta
				{ thePtrMedia->nuGradient(currLoc + .5*dL*tanNext, .125*dL)
				- thePtrMedia->nuGradient(currLoc - .5*dL*tanPrev, .125*dL)
				};
			BiVector const bivN{ nuDelta * wedge(tVal, (uVal)) };
			BiVector const bivT{ nuVal * wedge(tanDelta, (uVal)) };
			BiVector const bivU{ nuVal * wedge(tVal, normDelta) };
			BiVector const bivNT{ (nuDelta) * wedge(tanDelta, (uVal)) };
			BiVector const bivNU{ (nuDelta) * wedge(tVal, normDelta) };
			BiVector const bivTU{ nuVal * wedge(tanDelta, normDelta) };
			BiVector const bivNTU{ (nuDelta) * wedge(tanDelta, normDelta) };
			value =
				- bivN - bivT - bivU
				- bivNT - bivNU - bivTU
				- bivNTU
				;
std::cout << "  bivN: " << io::fixed(bivN) << '\n';
std::cout << "  bivT: " << io::fixed(bivT) << '\n';
std::cout << "  bivU: " << io::fixed(bivU) << '\n';
std::cout << " bivNT: " << io::fixed(bivNT) << '\n';
std::cout << " bivNU: " << io::fixed(bivNU) << '\n';
std::cout << " bivTU: " << io::fixed(bivTU) << '\n';
std::cout << "bivNTU: " << io::fixed(bivNTU) << '\n';

			// Estimate finite difference based on derivative values
			double const dnuDelta{ dL * (gDot*tVal).theSca[0] };
			Vector const dtanDelta{ dL * tDot };
			Vector const dnormDelta{ dL * uDot };
std::cout << "  dnuDelta: " << io::fixed(dnuDelta) << '\n';
std::cout << " dtanDelta: " << io::fixed(dtanDelta) << '\n';
std::cout << "dnormDelta: " << io::fixed(dnormDelta) << '\n';

std::cout << "        dL: " << io::fixed(dL) << '\n';
std::cout << "   nuDelta: " << io::fixed(nuDelta) << '\n';
std::cout << "  tanDelta: " << io::fixed(tanDelta) << '\n';
std::cout << " normDelta: " << io::fixed(normDelta) << '\n';


			return  value;
		}

	}; // DiffEq

	//! Check numerical path against candidate differential equation
	void
	test0
		( std::ostringstream & oss
		)
	{
		// Media with Gaussian blobs in it
		std::shared_ptr<tst::MediaBlobs> const ptrMedia
			{ std::make_shared<tst::MediaBlobs>() };

		// generate path by numerical propagation (essentially Euler's method)
		constexpr std::size_t reserveNodeSize{ 8u };
		std::shared_ptr<aply::ray::Path>
			const ptrPath{ numericalPath(ptrMedia, reserveNodeSize) };

		// evaluate equation at each node
oss << "Failure: Need to resolve what is expected from differential eqn\n";
		DiffEq const equation{ ptrMedia };
		std::vector<aply::ray::Node> const & nodes = ptrPath->theNodes;
		std::size_t const numNodes{ nodes.size() };
		if (3u < numNodes)
		{
			bool isFirst{ true };
			for (std::size_t nn{0u} ; nn < numNodes ; ++nn)
			{
std::cout << '\n';
				aply::ray::Node const & node = nodes[nn];
				std::cout << node.infoBrief() << std::endl;

				using namespace engabra::g3;
				BiVector const gap{ equation(node) };
// TODO adjust tolerance
				if (! (magnitude(gap) < 1.e-6))
				{
					if (! isFirst)
					{
					//	break;
					}
					isFirst = false;

					if (isFirst)
					{
						oss << "Failure of gap test\n";
						oss << "first gap: " << io::fixed(gap, 2u, 9u) << '\n';
					}
				}
			}
		}
	}
}


/*! \brief Compare numerical path with differential equation evaluation.
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);

	return tst::finish(oss);
}

