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
#include <vector>


namespace
{
	using namespace engabra::g3;

	constexpr double sNuEarth{ 1.010 }; // big (actual near 1.000273 STP)
	constexpr double sNuSpace{ 1.000 };
	constexpr double sRadEarth{ 6370.e3 };
	constexpr double sRadSpace{ sRadEarth + 100.e3 };

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


		double const theAlpha{};
		double const theBeta{};


		//! Construct model to match environment constants
		AtmModel
			()
			: theAlpha{ alpha(sNuEarth, sNuSpace, sRadEarth, sRadSpace) }
			, theBeta{ beta(sNuEarth, sNuSpace, sRadEarth, sRadSpace) }
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
		, double const & rBeg = sRadEarth
		, double const & rEnd = sRadSpace
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

	}; // Node

	//! Ray propagation functions
	struct Propagator
	{
		AtmModel const theAtm{};
		double const theStepDist{ null<double>() };

		//! True if this instance is valid
		inline
		bool
		isValid
			() const
		{
			return engabra::g3::isValid(theStepDist);
		}

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

		/*! Perform forward integration step by step.
		 *
		 * Essentially is Euler's method for integration of the ray path
		 * (with all attendent pitfalls).
		 *
		 */
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
				nodes.reserve(nomSize);

				// start with initial conditions
				Vector tPrev{ tBeg };
				Vector rCurr{ rBeg };
				double nuPrev{ null<double>() };

				// propagate until path approximate reaches requested length
				double length{ 0. };
				while (length < nominalLength)
				{
					// propagate ray through next step
					double nuNext{ null<double>() };
					Vector const tNext{ nextTangent(tPrev, rCurr, &nuNext) };
					Vector const rNext{ nextLocation(rCurr, tNext) };

					// record information for this node
					Node const node
						{ tPrev
						, nuPrev
						, rCurr
						, nuNext
						, tNext
						};
					nodes.emplace_back(node);

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
			<< " " << std::setw(3) << ndx
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
	AtmModel const atm;

	std::cout << "nu(sRadEarth): " << atm.nuValue(sRadEarth*e3) << '\n';
	std::cout << "nu(sRadSpace): " << atm.nuValue(sRadSpace*e3) << '\n';
	std::cout << "    sRadEarth: " << io::fixed(sRadEarth) << '\n';
	std::cout << "    sRadSpace: " << io::fixed(sRadSpace) << '\n';

	double const delta{ 10.e3 };
	double const nominalLength{ magnitude(sRadSpace - sRadEarth) };

	std::cout << "delta: " << io::fixed(delta) << '\n';

	// initial conditions
	Vector const tBeg{ direction(e1 + e3) };
	Vector const rBeg{ sRadEarth * e3 };

	// propagate ray forward until maxLength
	Propagator const prop{ atm, delta };
	std::vector<Node> const nodes{ prop.nodePath(tBeg, rBeg, nominalLength) };

	// report results
	std::size_t ndx{ 0u };
	for (Node const & node : nodes)
	{
		report(std::cout, node, ndx++);
	}

}

