//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


//! \file Atmospheric Refraction Example.

/*! \brief Demonstrate refraction path trace for exponential atmosphere.
 *
 * This is an example used to drive development of various software
 * classes in this project.
 *
 */


#include <Engabra>
#include <g3opsMul.hpp>

#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>


namespace
{
	using namespace engabra::g3;

	//! Data parameters for Planet atmosphere
	struct Planet
	{
		double const theNuGround{ null<double>() };
		double const theNuSpace{ null<double>() };
		double const theRadGround{ null<double>() };
		double const theRadSpace{ null<double>() };

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
			oss << " theNuGround: " << io::fixed(theNuGround);
			oss << '\n';
			oss << "  theNuSpace: " << io::fixed(theNuSpace);
			oss << '\n';
			oss << "theRadGround: " << io::fixed(theRadGround);
			oss << '\n';
			oss << " theRadSpace: " << io::fixed(theRadSpace);
			return oss.str();
		}


	}; // Planet

	//! Parameters for Earth atmosphere
	static Planet const sEarth
		{
		  1.000273 // nuGround // (for earth, actual is 1.000273 at STP)
		, 1.000    // nuSpace
		, 6370.e3  // radiusGround
		, 6470.e3  // radiusSpace
		};

	//! Atmospheric model : nu = alpha*exp(-beta*radius)
	struct AtmModel
	{
		//! Classic exponential decay model
		inline
		static
		double
		expValue
			( double const & amp
			, double const & scale
			, double const & radius
			)
		{
			return amp * std::exp(-scale * radius);
		}

		//! Exponential decay constant
		inline
		static
		double
		beta
			( double const & v1
			, double const & v2
			, double const & r1
			, double const & r2
			)
		{
			return { log(v1 / v2) / (r2 - r1) };
		}

		//! Exponential decay magnitude
		inline
		static
		double
		alpha
			( double const & v1
			, double const & v2
			, double const & r1
			, double const & r2
			)
		{
			double const bNeg{ - beta(v1, v2, r1, r2) };
			return
				{ (v1 + v2) / (std::exp(bNeg * r1) + std::exp(bNeg * r2)) };
		}


		double const theAlpha{ null<double>() };
		double const theBeta{ null<double>() };

		//! Construct an invalid instance
		inline
		AtmModel
			()
			: theAlpha{ null<double>() }
			, theBeta{ null<double>() }
		{
		}

		//! Construct model to match environment constants
		inline
		AtmModel
			( Planet const & planet
			)
			: theAlpha
				{ alpha
					( planet.theNuGround
					, planet.theNuSpace
					, planet.theRadGround
					, planet.theRadSpace
					)
				}
			, theBeta
				{ beta
					( planet.theNuGround
					, planet.theNuSpace
					, planet.theRadGround
					, planet.theRadSpace
					)
				}
		{ }

		//! Index of refraction value at vector location rVec
		inline
		double
		nuValue
			( Vector const & rVec
			) const
		{
			double const rMag{ magnitude(rVec) };
			return expValue(theAlpha, theBeta, rMag);
		}

		//! Gradient (approximation) of Index of refraction at location rVec
		inline
		Vector
		nuGradient
			( Vector const & rVec
			) const
		{
			double const mag{ magnitude(rVec) };
			double const del
				{ mag * std::sqrt(std::numeric_limits<double>::epsilon()) };
			double const dNu1
				{ nuValue(rVec + del * e1) - nuValue(rVec - del * e1) };
			double const dNu2
				{ nuValue(rVec + del * e2) - nuValue(rVec - del * e2) };
			double const dNu3
				{ nuValue(rVec + del * e3) - nuValue(rVec - del * e3) };
			double const scale{ 1. / (2. * del) };
			return { scale * Vector{ dNu1, dNu2, dNu3 } };
		}

		//! Unitary direction assocsiated with atm.nuGradient
		inline
		Vector
		gradDir
			( Vector const & rVec
			) const
		{
			return direction(nuGradient(rVec));
		}

	}; // AtmModel

	//! Sampling of atm.nuValue() values
	inline
	std::vector<double>
	nuProfile
		( double const & delta
		, AtmModel const & atm
		, double const & rBeg = sEarth.theRadGround
		, double const & rEnd = sEarth.theRadSpace
		)
	{
		std::vector<double> nus;
		// nus.reserve();
		for (double rad{rBeg} ; rad < rEnd ; rad += delta)
		{
			double const nu{ atm.nuValue(rad * e3) };
			nus.emplace_back(nu);
		}
		return nus;
	}

	//! Update tangent direction based on current step info
	inline
	Vector
	refractedTangent
		( Vector const & tPrev
		, Vector const & rCurr
		, double const & nuPrev
		, double const & nuNext
		, AtmModel const & atm
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
			Vector const uCurr{ atm.gradDir(rCurr) };

			// compute refraction bivector at current location
			BiVector const bivCurr{ (nuPrev/nuNext) * (tPrev*uCurr).theBiv };

			// Use radicand value to select between propagation methods
			// (total internal reflection or refraction)
			double const radicand{ 1. -magSq(bivCurr) };
			if (radicand < 0.)
			{
				// (internal) reflection case
				tNext = (uCurr * tPrev * uCurr).theVec;
			}
			else
			{
				// refraction case
				if (nuPrev < nuNext)
				{
					Spinor const spinU{  std::sqrt(radicand), bivCurr };
					tNext = (spinU * uCurr).theVec;
				}
				else if (nuNext < nuPrev)
				{
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

		return tNext;
	}

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

	//! Ray propagation functions
	struct Propagator
	{
		AtmModel const theAtm{};
		double const theStepDist{ null<double>() };

	private:

		//! Estimate next tangent based on local atmospheric refraction
		inline
		Vector
		nextTangent
			( Vector const & tPrev
			, Vector const & rCurr
			, double * const & ptrNuNext = nullptr
			) const
		{
			Vector tNext{ tPrev }; // implicit - needs iteration

			// incomming state is fixed
			double const nuPrev
				{ theAtm.nuValue(rCurr - .5*theStepDist*tPrev) };

			// iterate on exiting path
			double nuNext{ null<double>() }; // return value if requested
			double difSq{ 2. };
			static double const tol{ std::numeric_limits<double>::epsilon() };
			while (tol < difSq)
			{
				// update refraction index to midpoint of predicted next
				// interval (along evolving next tangent direction).
				nuNext = theAtm.nuValue(rCurr + .5*theStepDist*tNext);
				Vector const tTemp
					{ refractedTangent(tPrev, rCurr, nuPrev, nuNext, theAtm) };

				// evaluate convergence of tangent direction
				difSq = magSq(tTemp - tNext);

				// candidate return value (if convergence test passes)
				tNext = tTemp;
			}

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

				// start with initial conditions
				Vector tPrev{ tBeg };
				Vector rCurr{ rBeg };
				double nuPrev{ null<double>() };

				// propagate until path approximate reaches requested length
				double length{ 0. };
				std::size_t loopNum{ 0u };
				while (length < nominalLength)
				{
					// propagate ray through next step
					double nuNext{ null<double>() };
					Vector const tNext{ nextTangent(tPrev, rCurr, &nuNext) };
					Vector const rNext{ nextLocation(rCurr, tNext) };

					if (nodes.empty() || (0u == (loopNum++ % stride)))
					{
						// record information for this node
						Node const node
							{ tPrev
							, nuPrev
							, rCurr
							, nuNext
							, tNext
							};
						nodes.emplace_back(node);
					}

					// update state for next node
					tPrev = tNext;
					rCurr = rNext;
					nuPrev = nuNext;

					// track length (to terminate propagation)
					length += theStepDist;
				}
			}

			return nodes;
		}

	}; // Propagator

	//! Put current position and tangent values to stream
	void
	inline
	report
		( std::ostream & strm
		, Vector const & tPrev
		, Vector const & rCurr
		, Vector const & tNext
		, std::size_t const & ndx
		)
	{
		strm
			<< " ndx: " << std::setw(9u) << ndx
			<< " " << "tPrev: " << io::fixed(tPrev, 2u)
			<< " " << "rCurr: " << io::fixed(rCurr, 8u)
			<< " " << "tNext: " << io::fixed(tNext, 2u)
			<< '\n';
	}

	//! Put current node data values to stream
	void
	inline
	report
		( std::ostream & strm
		, Node const & node
		, std::size_t const & ndx
		)
	{
		report(strm, node.thePrevTan, node.theCurrLoc, node.theNextTan, ndx);
	}


} // [anon]

int
main
	()
{
	AtmModel const atm(sEarth);

	double const delta{ 10.e3 };
	double const nominalLength{ 150.e3 };
	//	{ magnitude(sEarth.theRadSpace - sEarth.theRadGround) };

	// initial conditions
	Vector const tBeg{ direction(e1 + e3) };
	Vector const rBeg{ sEarth.theRadGround * e3 };

	for (double delta{100000.} ; .000001 < delta ; delta /=10.)
	{
		Propagator const prop{ atm, delta };
		std::vector<Node> const fwdNodes
			{ prop.nodePath(tBeg, rBeg, nominalLength) };
		std::cout << " delta: " << io::fixed(delta, 7u, 6u) << " ";
		report(std::cout, fwdNodes.back(), fwdNodes.size());
	}

return 0;

	// propagate ray forward until maxLength
	Propagator const prop{ atm, delta };
	std::vector<Node> const fwdNodes
		{ prop.nodePath(tBeg, rBeg, nominalLength) };

	// report results
	std::cout << sEarth.infoString("sEarth") << '\n';
	for (std::size_t nn{0u} ; nn < fwdNodes.size() ; ++nn)
	{
		report(std::cout, fwdNodes[nn], nn);
	}
	/*
	std::cout << "        delta: " << io::fixed(delta) << '\n';
	std::cout << "nominalLength: " << io::fixed(nominalLength) << '\n';
	report(std::cout, fwdNodes.back(), fwdNodes.size());
	*/

}

/* -- TODO good unit test - checking reversal

	Node const & endNode = fwdNodes.back();
	Node const revNode{ endNode.reversed() };

	std::cout << " end node: \n" << endNode.infoString() << '\n';
	std::cout << " rev node: \n" << revNode.infoString() << '\n';

	// propagate backward
	std::vector<Node> const revNodes
		{ prop.nodePath
			(revNode.thePrevTan, revNode.theCurrLoc, nominalLength)
		};

	for (std::size_t nn{0u} ; nn < revNodes.size() ; ++nn)
	{
		report(std::cout, revNodes[nn], nn);
	}
*/
