//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//

#ifndef Refraction_ray_INCL_
#define Refraction_ray_INCL_

//! \file Ray propagation simulation functions.


#include "env.hpp"

#include <Engabra>

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>


/*! \brief Functions and classes for simulation of ray propagation
 */
namespace ray
{
	using namespace engabra::g3;

	//! Data relevant to an individual node
	struct Node
	{
		Vector const thePrevTan;
		double const thePrevNu;
		Vector const theCurrLoc;
		double const theNextNu;
		Vector const theNextTan;

		//! Descriptive information about this instance
		inline
		std::string
		infoBrief
			( std::string const & title = {}
			) const
		{
			std::ostringstream oss;
			if (! title.empty())
			{
				oss << title << '\n';
			}
			oss << " tan ";
			oss << io::fixed(thePrevTan, 3u, 6u);
			oss << " nu ";
			oss << io::fixed(thePrevNu, 3u, 6u);
			oss << " loc ";
			oss << io::fixed(theCurrLoc, 3u, 6u);
			oss << " nu ";
			oss << io::fixed(theNextNu, 3u, 6u);
			oss << " tan ";
			oss << io::fixed(theNextTan, 3u, 6u);
			return oss.str();
		}

		//! Descriptive information about this instance
		inline
		std::string
		infoString
			( std::string const & title = {}
			) const
		{
			std::ostringstream oss;
			if (! title.empty())
			{
				oss << title << '\n';
			}
			oss << "thePrevTan: " << io::fixed(thePrevTan, 8u, 6u);
			oss << '\n';
			oss << " thePrevNu: " << io::fixed(thePrevNu, 8u, 6u);
			oss << '\n';
			oss << "theCurrLoc: " << io::fixed(theCurrLoc, 8u, 6u);
			oss << '\n';
			oss << " theNextNu: " << io::fixed(theNextNu, 8u, 6u);
			oss << '\n';
			oss << "theNextTan: " << io::fixed(theNextTan, 8u, 6u);
			return oss.str();
		}

		//! Node associated with reversing direction of evolution
		inline
		Node
		reversed
			() const
		{
			return Node
				{ -theNextTan
				,  theNextNu
				,  theCurrLoc
				,  thePrevNu
				, -thePrevTan
				};
		}

	}; // Node

	//! Update tangent direction based on current step info
	inline
	Vector
	refractedTangent
		( Vector const & tPrev //!< Must be unit length
		, Vector const & rCurr
		, double const & nuPrev
		, double const & nuNext
		, env::IndexVolume const * const & ptObj
		)
	{
		Vector tNext{ null<Vector>() };

		// check validity
		// if (isValid(rCurr) && isValid(tPrev))
		// ... just assume valid for now
		{
			// Default corresponds with case (nuPrev == nuNext)
			tNext = tPrev;

			// get index field gradient
			Vector const gCurr{ ptObj->nuGradient(rCurr) };
			double const gMagSq{ magSq(gCurr) };
//std::cout << "  gCurr: " << io::fixed(gCurr, 1u, 18u) << std::endl;
//std::cout << "  gMagSq: " << io::fixed(gMagSq, 1u, 18u) << std::endl;
			if (std::numeric_limits<double>::epsilon() < gMagSq)
			{
				// unit direction of gradient
				// TODO - probably can redo formulae to work with full grad?
				Vector const uCurr{ (1./std::sqrt(gMagSq)) * gCurr };

				// compute refraction bivector at current location
				BiVector const bivCurr
					{ (nuPrev/nuNext) * (tPrev*uCurr).theBiv };

				// Use radicand value to select between propagation methods
				// (total internal reflection or refraction)
				double const radicand{ 1. - magSq(bivCurr) };

/*
std::cout << "  nuP/N: " << io::fixed(nuPrev/nuNext, 1u, 18u) << std::endl;
std::cout << "  tPrev: " << io::fixed(tPrev, 1u, 18u) << std::endl;
std::cout << "  uCurr: " << io::fixed(uCurr, 1u, 18u) << std::endl;
std::cout << "  bivCurr: " << io::fixed(bivCurr, 1u, 18u) << std::endl;
std::cout << "  radicand: " << io::fixed(radicand, 1u, 18u) << std::endl;
*/

constexpr bool showIt{ false };
if (showIt)
{
std::cout
	<< "nuPrev: " << io::fixed(nuPrev)
	<< "  "
	<< "nuNext: " << io::fixed(nuNext)
	<< "  "
	;
if (nuPrev > nuNext)
	std::cout << "  ->  ";
if (nuPrev < nuNext)
	std::cout << "   <- ";
if (nuPrev == nuNext)
	std::cout << "   == ";
std::cout
	<< "  tPrev: " << io::fixed(tPrev, 1u, 18u)
	<< "  rCurr: " << io::fixed(rCurr, 1u, 18u)
	<< '\n';
std::cout << "  uCurr: " << io::fixed(uCurr, 1u, 18u) << std::endl;
std::cout << "  radicand: " << io::fixed(radicand, 1u, 18u) << std::endl;
}
				if (radicand < 0.)
				{
					// (internal) reflection case
					tNext = -(uCurr * tPrev * uCurr).theVec;
std::cout << "\n###### REFLECTION #########\n";
				}
				else
				{
					// refraction case
					if (nuPrev < nuNext)
					{
//std::cout << "Spin-A\n";
						Spinor const spinU{  std::sqrt(radicand), bivCurr };
						tNext = (spinU * uCurr).theVec;
					}
					else if (nuNext < nuPrev)
					{
//std::cout << "Spin-B\n";
						Spinor const spinU{ -std::sqrt(radicand), bivCurr };
						tNext = (spinU * uCurr).theVec;
					}
				//	// Default case
				//	else // nuPrev == nuNext
				//	{
				//		tNext = tPrev; // could initialize with this?
				//	}
				}
			}
		}

		return tNext;
	}

	//! Ray propagation functions
	struct Propagator
	{
		env::IndexVolume const * const thePtObj{ nullptr };
		double const theStepDist{ null<double>() };

	private:

		//! Estimate next tangent based on local object refraction
		inline
		Vector
		nextTangent
			( Vector const & tPrev //!< Must be unit length
			, Vector const & rCurr
			, double * const & ptrNuNext = nullptr
			) const
		{
			Vector tNext{ tPrev }; // implicit - needs iteration

			// incomming state is fixed
			double const nuPrev
				{ thePtObj->nuValue(rCurr - .5*theStepDist*tPrev) };

			// iterate on exiting path
			double nuNext{ null<double>() }; // return value if requested
			double difSq{ 2. };
			static double const tol{ std::numeric_limits<double>::epsilon() };
			std::size_t numLoop{ 0u };
constexpr std::size_t maxLoop{ 10u };
//TODO does this converge?
//std::cout << "---- looping\n";
			while ((tol < difSq) && (numLoop++ < maxLoop))
			{
				// update refraction index to midpoint of predicted next
				// interval (along evolving next tangent direction).
				nuNext = thePtObj->nuValue(rCurr + .5*theStepDist*tNext);
				Vector const tTemp
					{ refractedTangent
						(tPrev, rCurr, nuPrev, nuNext, thePtObj)
					};

/*
Vector const rPrev{ rCurr - .5*theStepDist*tPrev };
Vector const rNext{ rCurr + .5*theStepDist*tNext };
std::cout
	<< "  rPrev: " << io::fixed(rPrev, 1u, 9u)
	<< "  rCurr: " << io::fixed(rCurr, 1u, 9u)
	<< "  rNext: " << io::fixed(rNext, 1u, 9u)
	<< "  tTemp: " << io::fixed(tTemp, 1u, 9u)
	<< std::endl;
*/

				// evaluate convergence of tangent direction
				difSq = magSq(tTemp - tNext);

				// candidate return value (if convergence test passes)
				tNext = tTemp;
			}
//std::cout << "--\n";

			// set for use in consumer code (if requested)
			if (ptrNuNext)
			{
				*ptrNuNext = nuNext;
			}

			return tNext;
		}

		//! Predicted next location stepsize units along tangent from rVec
		inline
		Vector
		nextLocation
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
		isValid
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
		traceNodes
			( Vector const & tBeg
			, Vector const & rBeg
			, Consumer * const & ptConsumer
			) const
		{
			if (isValid())
			{
				// start with initial conditions
				Vector tPrev{ direction(tBeg) };
				Vector rCurr{ rBeg };
				double nuPrev{ null<double>() };

				// propagate until path approximate reaches requested length
				while (ptConsumer->size() < ptConsumer->capacity())
				{
					// propagate ray through next step
					double nuNext{ null<double>() };
					Vector const tNext{ nextTangent(tPrev, rCurr, &nuNext) };
					Vector const rNext{ nextLocation(rCurr, tNext) };

					// Provide intermediate info to consumer
					ptConsumer->emplace_back
						(Node{ tPrev, nuPrev, rCurr, nuNext, tNext });

					// update state for next node
					tPrev = tNext;
					rCurr = rNext;
					nuPrev = nuNext;
				}
std::cout << "\n####\n#### done tracing\n####\n";
			}
		}

/*TODO*/
		//! \brief Return nodes in a std::vector structure.
		inline
		std::vector<Node>
		nodePath
			( Vector const & tBeg
			, Vector const & rBeg
			, double const & nominalLength
			) const
		{
			std::vector<Node> nodes;

			if (isValid())
			{
				// estimate return structure size and allocate space
				constexpr double padFactor{ 9./8. }; // about 12% extra
				double const dubSize{ padFactor * nominalLength / theStepDist };
				std::size_t const nomSize{ static_cast<std::size_t>(dubSize) };

				std::size_t stride{ 1u };
				std::size_t saveSize{ nomSize };
				if (1000u < nomSize)
				{
					stride = nomSize / 1000u;
				}
				std::size_t const useSize{ nomSize / stride + 10u };
				nodes.reserve(useSize);

				traceNodes(tBeg, rBeg, &nodes);
std::cout << "nodes.size; " << nodes.size() << std::endl;
			}

			return nodes;
		}
/*TODO*/

	}; // Propagator

	//! Put current position and tangent values to stream
	std::string
	inline
	nodeStateInfo
		( Vector const & tPrev
		, Vector const & rCurr
		, Vector const & tNext
		, std::size_t const & ndx
		)
	{
		std::ostringstream oss;
		oss
			<< " ndx: " << std::setw(9u) << ndx
			<< " " << "tPrev: " << io::fixed(tPrev, 2u)
			<< " " << "rCurr: " << io::fixed(rCurr, 8u)
			<< " " << "tNext: " << io::fixed(tNext, 2u)
			;
		return oss.str();
	}

	//! Put current node data values to stream
	std::string
	inline
	nodeInfo
		( ray::Node const & node
		, std::size_t const & ndx
		)
	{
		return nodeStateInfo
			(node.thePrevTan, node.theCurrLoc, node.theNextTan, ndx);
	}



} // [ray]

#endif // Refraction_ray_INCL_

