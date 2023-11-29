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
	snelsTangent
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

	//! Estimate next tangent direction based on local atmospheric refraction
	inline
	Vector
	nextTangent
		( Vector const & tPrev
		, Vector const & rCurr
		, AtmModel const & atm
		, double const & delta
		)
	{
		Vector tNext{ tPrev }; // implicit - needs iteration

		// incomming state is fixed
		double const nuPrev{ atm.nuValue(rCurr - .5*delta*tPrev) };

		// iterate on exiting path
		double difSq{ 2. };
		static double const tol{ std::numeric_limits<double>::epsilon() };
		while (tol < difSq)
		{
			double const nuNext{ atm.nuValue(rCurr + .5*delta*tNext) };
			Vector const tTemp
				{ snelsTangent(tPrev, rCurr, nuPrev, nuNext, atm) };
			difSq = magSq(tTemp - tNext);
			tNext = tTemp;
		}

		return tNext;
	}

	//! Predicted next location delta units along tangent from rVec
	inline
	Vector
	nextLocation
		( Vector const & rVec
		, Vector const & tVec
		, double const & delta
		)
	{
		return { rVec + delta * tVec };
	}

	//! Put current position and tangent values to stream
	void
	inline
	update
		( std::ostream & strm
		, Vector const & tVec
		, Vector const & rVec
		, std::size_t const & ndx
		)
	{
		strm
			<< " " << std::setw(3) << ndx
			<< " " << "tan: " << io::fixed(tVec, 2u)
			<< " " << "loc: " << io::fixed(rVec, 8u)
			<< '\n';
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

	Vector const tBeg{ direction(e1 + e3) };
	Vector const rBeg{ sRadEarth * e3 };
update(std::cout, tBeg, rBeg, 0u);

	Vector const t1{ nextTangent(tBeg, rBeg, atm, delta) };
	Vector const r1{ nextLocation(rBeg, t1, delta) };
update(std::cout, t1, r1, 1u);

	Vector const t2{ nextTangent(t1, r1, atm, delta) };
	Vector const r2{ nextLocation(r1, t2, delta) };
update(std::cout, t2, r2, 2u);

}

