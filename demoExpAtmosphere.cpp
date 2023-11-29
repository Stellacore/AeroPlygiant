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
		{
		}

		//! Value of index at scalar radial distance
		inline
		double
		nuValue
			( double const & radMag
			) const
		{
			return expValue(theAlpha, theBeta, radMag);
		}

		//! Index of refraction value at vector location rVec
		inline
		double
		nuValue
			( Vector const & rVec
			) const
		{
			double const rMag{ magnitude(rVec) };
			return nuValue(rMag);
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
	nextTangent
		( Vector const & rCurr
		, Vector const & tCurr
		, double const & nuCurr
		, double const & nuNext
		, AtmModel const & atm
		)
	{
		Vector tNext{ null<Vector>() };

		// check validity
		// if (isValid(rCurr) && isValid(tCurr))
		// ... just assume valid for now
		{
			// Default corresponds with case (nuCurr == nuNext)
			tNext = tCurr;

			// get index field gradient
			Vector const uCurr{ atm.gradDir(rCurr) };

			// compute refraction bivector at current location
			BiVector const bivCurr{ (nuCurr/nuNext) * (tCurr*uCurr).theBiv };

			// Use radicand value to select between propagation methods
			// (total internal reflection or refraction)
			double const radicand{ 1. -magSq(bivCurr) };
			if (radicand < 0.)
			{
				// (internal) reflection case
				tNext = (uCurr * tCurr * uCurr).theVec;
			}
			else
			{
				// refraction case
				if (nuCurr < nuNext)
				{
					Spinor const spinU{  std::sqrt(radicand), bivCurr };
					tNext = (spinU * uCurr).theVec;
				}
				else if (nuNext < nuCurr)
				{
					Spinor const spinU{ -std::sqrt(radicand), bivCurr };
					tNext = (spinU * uCurr).theVec;
				}
			//	// Default case
			//	else // nuCurr == nuNext
			//	{
			//		tNext = tCurr; // could initialize with this?
			//	}
			}
		}

		return tNext;
	}

	//! Representation of location and tangent at point along path
	struct NodeRT
	{
		Vector const theLoc{ null<Vector>() };
		Vector const theTan{ null<Vector>() };

		//! Location of point on ray
		inline
		Vector const &
		loc
			() const
		{
			return theLoc;
		}

		//! Path tangent (unit) direction at loc()
		inline
		Vector const &
		tan
			() const
		{
			return theTan;
		}

	};

	//! Next node (location and exiting tangent : <rVec,tVec>)
	inline
	NodeRT
	nextNode
		( NodeRT const & rtCurr
		, double const & nuCurr
		, double const & nuNext
		, AtmModel const & atm
		, double const & delta
		)
	{
		Vector const & rCurr = rtCurr.loc();
		Vector const & tCurr = rtCurr.tan();
		Vector const tNext{ nextTangent(rCurr, tCurr, nuCurr, nuNext, atm) };
		Vector const rNext{ rCurr + delta * tNext };
		return { rNext, tNext };
	}

	//! Put current position and tangent values to stream
	void
	inline
	update
		( std::ostream & strm
		, Vector const & rVec
		, Vector const & tVec
		, std::size_t const & ndx
		)
	{
		strm
			<< " " << std::setw(3) << ndx
			<< " " << "loc: " << io::fixed(rVec, 8u)
			<< " " << "tan: " << io::fixed(tVec, 2u)
			<< '\n';
	}

	//! Put current (position, tangent) pari values to stream
	void
	inline
	update
		( std::ostream & strm
		, NodeRT const & rtPair
		, std::size_t const & ndx
		)
	{
		update(strm, rtPair.loc(), rtPair.tan(), ndx);
	}

} // [anon]


int
main
	()
{
	AtmModel const atm;

	std::cout << "nu(sRadEarth): " << atm.nuValue(sRadEarth) << '\n';
	std::cout << "nu(sRadSpace): " << atm.nuValue(sRadSpace) << '\n';
	std::cout << "    sRadEarth: " << io::fixed(sRadEarth) << '\n';
	std::cout << "    sRadSpace: " << io::fixed(sRadSpace) << '\n';

	double const delta{ 10.e3 };
	std::cout << "delta: " << io::fixed(delta) << '\n';

	Vector const r0{ sRadEarth * e3 };
	Vector const t0{ direction(e1 + e3) };
	NodeRT const rt0{ r0, t0 };
update(std::cout, rt0, 0u);

	double const nu0{ atm.nuValue(r0 - .5*delta*t0) };
	double const nu1{ atm.nuValue(r0 + .5*delta*t0) };
	NodeRT const rt1{ nextNode(rt0, nu0, nu1, atm, delta) };
update(std::cout, rt1, 1u);

	double const nu2{ atm.nuValue(rt1.loc() + .5*delta*rt1.tan()) };
	NodeRT const rt2{ nextNode(rt1, nu1, nu2, atm, delta) };
update(std::cout, rt2, 2u);

}

