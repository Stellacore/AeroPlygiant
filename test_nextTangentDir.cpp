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
 * \brief Unit test for class nextTangentDir
 *
 */


#include "tst.hpp"

#include "ray.hpp"


namespace
{
	using namespace engabra::g3;

	//! Rotate normDir in rayPlane by inAngleMag.
	Vector
	rotatedDir
		( Vector const & normDir
		, double const inAngleMag
		, BiVector const & rayPlane
		)
	{
		BiVector const angle{ inAngleMag * direction(rayPlane) };
		Spinor const spin{ exp(.5 * angle) };
		// active convention for rotation
		Vector const inDir{ (reverse(spin) * normDir * spin).theVec };
		return inDir;
	}

	//! Test case sample
	struct Sample
	{
		double const theNuIn; // on incoming side
		Vector const theNormDir; // from In to Out
		double const theNuOt;  // on outgoing side

		BiVector const theInPlane; // plane of incidence

		double const theInAngle;
		double const theOtAngle;

		explicit
		Sample
			( Vector const & normDir
			, Vector const & orthDir // norm/orthDir define ray plane
			, double const & nuIn
			, double const & nuOt
			, double const & inAngle
			, double const & otAngle
			)
			: theNuIn{ nuIn }
			, theNormDir{ direction(normDir) }
			, theNuOt{ nuOt }
			, theInPlane{ direction((normDir * orthDir).theBiv) }
			, theInAngle{ inAngle }
			, theOtAngle{ otAngle }
		{ }

		Vector
		tanIn
			() const
		{
			Vector const dir
				{ rotatedDir(theNormDir, theInAngle, theInPlane) };
			return dir;
		}

		Vector
		tanOt
			() const
		{
			Vector const dir
				{ rotatedDir(theNormDir, theOtAngle, theInPlane) };
			return dir;
		}

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
			oss << " theNuIn: " << theNuIn;
			oss << '\n';
			oss << "theNormDir: " << theNormDir;
			oss << '\n';
			oss << " theNuOt: " << theNuOt;
			oss << '\n';
			oss << "theInPlane: " << theInPlane;
			oss << '\n';
			oss << "theInAngle: " << theInAngle << std::endl;
			oss << '\n';
			oss << "theOtAngle: " << theOtAngle;
			oss << '\n';
			oss << " tanIn(): " << tanIn();
			oss << '\n';
			oss << " tanOt(): " << tanOt();
			return oss.str();
		};

	};

	//! Run test on refracted result
	void
	checkRefract
		( std::ostringstream & oss
		, std::pair<Vector, ray::DirChange> const & gotDirChange
		, std::pair<Vector, ray::DirChange> const & expDirChange
		, std::string const & tname
		)
	{
		// check tangent
		Vector const & expTanNext = expDirChange.first;
		Vector const & gotTanNext = gotDirChange.first;
		constexpr double tol{ 2. * std::numeric_limits<double>::epsilon() };
		if (! nearlyEquals(gotTanNext, expTanNext, tol))
		{
			Vector const difTanNext{ gotTanNext - expTanNext };
			oss << "Failure of next tangent '" << tname << "' test\n";
			oss << "expTanNext: " << expTanNext << '\n';
			oss << "gotTanNext: " << gotTanNext << '\n';
			oss << "difTanNext: " << io::fixed(difTanNext, 1u, 18u) << '\n';
		}

		// check change indicator
		ray::DirChange const & expChange = expDirChange.second;
		ray::DirChange const & gotChange = gotDirChange.second;
		if (! (expChange == gotChange))
		{
			oss << "Failure of refraction change '" << tname << "' test\n";
			std::cout << "expChange: " << ray::nameFor(expChange) << '\n';
			std::cout << "gotChange: " << ray::nameFor(gotChange) << '\n';
		}
	}

	// Config:
	//   Air
	//      e3 (normal)
	// Glass
	constexpr double sNuOt{ 1.  };
	constexpr double sNuIn{ 1.5 };
	Vector const sGradDir{ -e3 };
	constexpr double sDubNumRefract{ 4. };
	static Vector const sUp{ e3 }; // normal to interface
	static Vector const sOrth{ e1 }; // orthogonal to normal
	//
	// up to critical angle expect refraction forward and reverse
	static double const sInAngMagC{ std::asin(sNuOt / sNuIn) }; // critical ang
	static double const sDelAng{ sInAngMagC / sDubNumRefract };


	//! Check computation of boundary layer propagation (before critical angle)
	std::string
	test0_RefractGlassToAir
		()
	{
		std::ostringstream oss;

		// [DoxyExample00]
		// [DoxyExample00]

		for (double inA{0.} ; (! (sInAngMagC < inA)) ; inA += sDelAng)
		{
			// Snel's law
			double const otSinAng{ (sNuIn/sNuOt) * std::sin(inA) };
			double const otA{ std::asin(otSinAng) };
			Sample const sample(sUp, sOrth, sNuIn, sNuOt, inA, otA);

			// forward direction (glass to air) (refraction up to critical)
			{
				std::pair<Vector, ray::DirChange> const expDirChange
					{ sample.tanOt(), ray::Diverged };

				// get computed result
				std::pair<Vector, ray::DirChange> const gotDirChange
					{ ray::nextTangentDir
						(sample.tanIn(), sNuIn, sGradDir, sNuOt)
					};

				// check result
				checkRefract(oss, gotDirChange, expDirChange, "FWD-A");
			}

			// reverse direction (air to glass) (refraction as well)
			{
				std::pair<Vector, ray::DirChange> const expDirChange
					{ -sample.tanIn(), ray::Converged };

				// get computed result
				std::pair<Vector, ray::DirChange> const gotDirChange
					{ ray::nextTangentDir
						(-sample.tanOt(), sNuOt, sGradDir, sNuIn)
					};

				// check result
				checkRefract(oss, gotDirChange, expDirChange, "REV-A");
			}

			// append test configuration info if any test fails
			if (! oss.str().empty())
			{
				oss << sample.infoString("\nFor TestSample-A") << std::endl;
				break;
			}

		} // angles less than or equal to forward crtical

		return oss.str();
	}

	//! Check computation of boundary layer propagation (after critical angle)
	std::string
	test1_ReflectGlassToAir
		()
	{
		std::ostringstream oss;

		// [DoxyExample00]
		// [DoxyExample00]

		// after critical angle: reflection forward and refraction reverse
		constexpr double eps{ std::numeric_limits<double>::epsilon() };
		for (double inA{sInAngMagC+eps} ; (! (piHalf < inA)) ; inA += sDelAng)
		{
			// Reflection
			double const otA{ pi - inA };
			Sample const sample(sUp, sOrth, sNuIn, sNuOt, inA, otA);

			// forward direction (glass to air) reflection after critical angle.
			{
				std::pair<Vector, ray::DirChange> const expDirChange
					{ sample.tanOt(), ray::Reflected };

				// get computed result
				std::pair<Vector, ray::DirChange> const gotDirChange
					{ ray::nextTangentDir
						(sample.tanIn(), sNuIn, sGradDir, sNuOt)
					};

				// check result
				checkRefract(oss, gotDirChange, expDirChange, "FWD-B");
			}

			// append test configuration info if any test fails
			if (! oss.str().empty())
			{
				oss << sample.infoString("\nFor TestSample-B") << std::endl;
				break;
			}

		} // angles after critical to piHalf

		return oss.str();
	}

	//! Check computation of boundary layer propagation (before critical angle)
	std::string
	test2_RefractAirToGlass
		()
	{
		std::ostringstream oss;

		// [DoxyExample00]
		// [DoxyExample00]

		// reverse indices from above for this test
		constexpr double nuIn{ sNuOt };
		constexpr double nuOt{ sNuIn };

		for (double inA{0.} ; (! (piHalf < inA)) ; inA += sDelAng)
		{
			// Snel's law
			double const otSinAng{ (nuIn/nuOt) * std::sin(inA) };
			double const otA{ std::asin(otSinAng) };
			// Note reversal of the test-sample normal direction
			Vector const down{ -sUp };
			Sample const sample(down, sOrth, nuIn, nuOt, inA, otA);

			/*
			std::cout << std::endl;
			std::cout << "  inA: " << inA << "  nuIn: " << nuIn << '\n';
			std::cout << "inDir: " << sample.tanIn() << '\n';
			std::cout << " norm: " << sample.theNormDir << '\n';
			std::cout << "  sUp: " << sUp << '\n';
			std::cout << "sGrad: " << sGradDir << '\n';
			std::cout << "otDir: " << sample.tanOt() << '\n';
			std::cout << "  otA: " << otA << "  nuOt: " << nuOt << '\n';
			*/

			// forward direction (glass to air) (refraction up to critical)
			{
				std::pair<Vector, ray::DirChange> const expDirChange
					{ sample.tanOt(), ray::Converged };

				// get computed result
				std::pair<Vector, ray::DirChange> const gotDirChange
					{ ray::nextTangentDir
						(sample.tanIn(), nuIn, sGradDir, nuOt)
					};

				// check result
				checkRefract(oss, gotDirChange, expDirChange, "FWD-C");
			}

			// append test configuration info if any test fails
			if (! oss.str().empty())
			{
				oss << sample.infoString("\nFor TestSample-C") << std::endl;
				break;
			}

		} // angles less than or equal to forward crtical

		return oss.str();
	}

}


/*! \brief Unit test for nextTangentDir
 */
int
main
	()
{
	std::ostringstream oss;

	oss << test0_RefractGlassToAir();
	oss << test1_ReflectGlassToAir();
	oss << test2_RefractAirToGlass();

	return tst::finish(oss);
}

