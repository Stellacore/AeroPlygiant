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
	constexpr double sPropStepDist{ 1./1024. };
	constexpr double sSaveStepDist{ 1./8. };

	//! Numerically generated path
	std::shared_ptr<aply::ray::Path>
	numericalPath
		( std::shared_ptr<tst::MediaBlobs> const & ptrMedia
		)
	{
		std::shared_ptr<aply::ray::Path> ptrPath{ nullptr };

		using namespace aply;
		using namespace engabra::g3;

		// Initial conditions
		ray::Start const start{ ray::Start::from(e1, zero<Vector>()) };

		// Allocate path
		ptrPath = std::make_shared<ray::Path>(start, sSaveStepDist);
		ptrPath->reserve(1128u);

		// Trace path through a 3D refraction environment
		ray::Propagator const prop{ ptrMedia.get(), sPropStepDist };
		prop.tracePath(ptrPath.get());

		return ptrPath;
	}

	struct DiffEq
	{
		std::shared_ptr<tst::MediaBlobs> const thePtrMedia{ nullptr };

		engabra::g3::MultiVector
		operator()
			( aply::ray::Node const & node
			) const
		{
			using namespace engabra::g3;

			// current location at which to evaluate diffeq
			Vector const & loc = node.theCurrLoc;

			// estimate point location tangent and index of refraction
			// values as average of before/after path sections
			Vector const tan{ .5 * (node.theNextTan + node.thePrevTan) };
			double const nu{ .5 * (node.theNextNu + node.thePrevNu) };

			// Estimate the derivative of the tangent (w.r.t. arc length)
			Vector const dTan
				{ (1./sPropStepDist) * (node.theNextTan - node.thePrevTan) };

			// Get the gradient vector at this location
			Vector const grad{ thePtrMedia->nuGradient(loc, sPropStepDist) };

			// Get the Hessian values at this location
			aply::math::Matrix const hess
				{ thePtrMedia->nuHessian(loc, sPropStepDist) };

			Vector const tgSum{ tan + grad };
			Vector const tgDif{ tan - grad };
			Vector const v1{ (1./nu) * (grad * tan).theSca[0] * grad };
			Vector const v2
				{ hess[0][0]*tan[0] + hess[0][1]*tan[1] + hess[0][2]*tan[2]
				, hess[1][0]*tan[0] + hess[1][1]*tan[1] + hess[1][2]*tan[2]
				, hess[2][0]*tan[0] + hess[2][1]*tan[1] + hess[2][2]*tan[2]
				};

			MultiVector const lhs1{ tgSum * dTan };
			MultiVector const lhs2{ dTan * tgDif };
			MultiVector const bigD{ 2. * (tan * (v1 + v2)).theBiv };

			MultiVector const lhs{ lhs1 + lhs2 };
			MultiVector const rhs{ bigD };

			MultiVector const eqn{ lhs - rhs };

std::cout << "lhs: " << io::fixed(lhs) << '\n';
std::cout << "rhs: " << io::fixed(rhs) << '\n';
//std::cout << "eqn: " << io::fixed(eqn) << '\n';

			return eqn;
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
		std::shared_ptr<aply::ray::Path>
			const ptrPath{ numericalPath(ptrMedia) };

		// evaluate equation at each node
		DiffEq const eqn{ ptrMedia };
		bool isFirst{ true };
		for (aply::ray::Node const & node : ptrPath->theNodes)
		{
std::cout << '\n';
			// Here, just display brief summary of individual node information
			std::cout << node.infoBrief() << std::endl;

			using namespace engabra::g3;
			MultiVector const gap{ eqn(node) };
std::cout << "gap: " << io::fixed(gap) << '\n';
			if (! (magnitude(gap) < 1.e-6))
			{
				if (! isFirst)
				{
					break;
				}
				isFirst = false;
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

oss << "Failure code me\n";
	return tst::finish(oss);
}

