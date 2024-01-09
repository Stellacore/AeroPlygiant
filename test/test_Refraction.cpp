//
// Original Code:
//  Copyright (c) 2007 Stellacore Corporation.
//  Donated to AeroPlygiant Open Source project 2024.
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
\brief Compare computed refraction with Manual of Photogrammetry table data.
*/


#include "tst.hpp"

#include "envPlanet.hpp"
#include "rayRefraction.hpp"

#include <Engabra>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <sstream>
#include <vector>


namespace
{

/* \brief Refraction reference data used for external validation.

Source: Manual of Photogrammetry, Forth Edition,
        Editor in Chief Chester C Slama
		American Society of Photogrammetry, Falls Church, Va
		1980, pg. 487

Refraction angle in uRad at 45-deg for ground high 0 km as function of
sensor high - source "3" (3rd Edition of M.o.P.)

H[km]    refraction[uRad]

   .5          4.9
  1.0          9.8
  1.5         14.9
  2.0         19.9
  2.5         25.0
  3.0         30.0
  3.5         35.0
  4.0         39.8
  4.5         44.6
  5.0         49.2
  5.5         53.6
  6.0         57.8
  6.5         61.9
  7.0         65.6
  7.5         69.2
  8.0         72.5
  8.5         75.5
  9.0         78.3
  9.5         80.8
 10.0         83.1

*/
namespace mop
{
	using AltDeviation = std::pair<double, double>;

	//!< Ref "namespace mop" for description and source.
	std::vector<AltDeviation> const altDevs
	{
		  {   .5e3 ,   4.9e-6 }
		, {  1.0e3 ,   9.8e-6 }
		, {  1.5e3 ,  14.9e-6 }
		, {  2.0e3 ,  19.9e-6 }
		, {  2.5e3 ,  25.0e-6 }
		, {  3.0e3 ,  30.0e-6 }
		, {  3.5e3 ,  35.0e-6 }
		, {  4.0e3 ,  39.8e-6 }
		, {  4.5e3 ,  44.6e-6 }
		, {  5.0e3 ,  49.2e-6 }
		, {  5.5e3 ,  53.6e-6 }
		, {  6.0e3 ,  57.8e-6 }
		, {  6.5e3 ,  61.9e-6 }
		, {  7.0e3 ,  65.6e-6 }
		, {  7.5e3 ,  69.2e-6 }
		, {  8.0e3 ,  72.5e-6 }
		, {  8.5e3 ,  75.5e-6 }
		, {  9.0e3 ,  78.3e-6 }
		, {  9.5e3 ,  80.8e-6 }
		, { 10.0e3 ,  83.1e-6 }
	};

} // [mop]


/*! \brief Check trivial cases
 *
 */
std::string
test0
	( std::ostringstream & oss
	)
{
	//
	// Test zero refraction for zero inclination angle
	//

	double const radEarth{ aply::env::sEarth.theRadGround };
	double const radGnd{ radEarth + 1000. };
	double const radSen{   radGnd + 1000. };

	aply::ray::Refraction const zeroRefract(0., radSen, radEarth);
	double const gotZero(zeroRefract.thetaAngleAt(radGnd));
	double const expZero(0.0);

	if (! engabra::g3::nearlyEquals(gotZero, expZero))
	{
		using engabra::g3::io::fixed;
		oss << "Failure of zero inclination angle test:" << std::endl;
		oss << "got: " << fixed(gotZero, 1u, 9u) << std::endl;
		oss << "exp: " << fixed(expZero, 1u, 9u) << std::endl;
	}

	return oss.str();
}

/*! \brief Test consistency of refraction deviation angles.
 */
void
checkRefDevAngles
	( std::ostream & oss
	, double const & gotRefDevAngle
	, double const & expRefDevAngle
	, double const & highSensor
	, double const & highGround
	, double const & tolAngleAbsolute
	)
{
	// Note 'absolute' version of comparisionf
	// (since relative errors on this data are enormous (couple percent)).
	using engabra::g3::nearlyEqualsAbs;
	if (! nearlyEqualsAbs(gotRefDevAngle, expRefDevAngle, tolAngleAbsolute))
	{
		using engabra::g3::io::fixed;
		double const difRefDevAngle{ gotRefDevAngle - expRefDevAngle };
		oss << "Failure of forward refraction angle test\n";
		oss << "    highSensor: " << fixed(highSensor, 5u, 3u)
			<< "  [m]\n";
		oss << "    highGround: " << fixed(highGround, 5u, 3u)
			<< "  [m]\n";
		oss << "expRefDevAngle: " << fixed(expRefDevAngle, 1u, 6u)
			<< "  From MoP {3rd Ed., pg487}\n";
		oss << "gotRefDevAngle: " << fixed(gotRefDevAngle, 1u, 6u)
			<< "  Using COESA1976 Atmosphere model\n";
		oss << "difRefDevAngle: " << fixed(difRefDevAngle, 1u, 6u) << '\n';
	}
}


/*! \brief Check integration of Gyer Eqn[12] for example high altitude use case.
 */
std::string
test1
	( std::ostringstream & oss
	)
{
	// Use for outputting errors
	oss.precision(20);

	// for convenience
	using engabra::g3::nearlyEquals;
	using engabra::g3::io::fixed;

	// [DoxyExample00]

	//
	// Example similar to mid oblique remote sensing geometry
	// Matches mop::AltDeviation validation data table described 
	// in namespace mop above which includes data from Manual
	// of Photogrammetry (3rd Ed, p 487).
	//
	constexpr double fwdLookAngle{ engabra::g3::piQtr }; // 45-deg off Nadir
	constexpr double highGround{     0. }; // [m] - "sea level" for MoP compare

	constexpr double highSensor{  9000. }; // [m] - a bit under 30k' (FL300)
	constexpr double expRefDevAngle{ .000078300 }; // from MoP table

	// convert to geocentric values for refraction computation
	double const radEarth{ aply::env::sEarth.theRadGround };
	double const radSen{ radEarth + highSensor };
	double const radGnd{ radEarth + highGround };

	aply::ray::Refraction const refract(fwdLookAngle, radSen, radEarth);
	double const thetaAtEnd{ refract.thetaAngleAt(radGnd) };

	double const gotRefDevAngle
		{ refract.angularDeviationFromStart(radGnd, thetaAtEnd) };

	// [DoxyExample00]

	// check consistency
	constexpr double tolAngle{ .000005 }; // about 1 arc second
	checkRefDevAngles
		( oss
		, gotRefDevAngle
		, expRefDevAngle
		, highSensor
		, highGround
		, tolAngle
		);

	return oss.str();
}

/*! \brief Check computations against MoP table of data.
 */
std::string
test2
	( std::ostringstream & oss
	)
{
	std::ostringstream rptResid;

	// if true, display residual to stdout
	constexpr bool showResiduals{ true };

	constexpr double fwdLookAngle{ engabra::g3::piQtr }; // 45-deg off Nadir
	constexpr double highGround{     0. }; // [m] - "sea level" for MoP compare

	// track differences for each table entry
	std::vector<double> resids{};
	resids.reserve(mop::altDevs.size());

	//! Loop over entire table
	for (mop::AltDeviation const & altDev : mop::altDevs)
	{
		double const & highSensor = altDev.first;
		double const & expRefDevAngle = altDev.second;

		// convert to geocentric values for refraction computation
		double const radEarth{ aply::env::sEarth.theRadGround };
		double const radSen{ radEarth + highSensor };
		double const radGnd{ radEarth + highGround };

		aply::ray::Refraction const refract(fwdLookAngle, radSen, radEarth);
		double const thetaAtEnd{ refract.thetaAngleAt(radGnd) };

		double const gotRefDevAngle
			{ refract.angularDeviationFromStart(radGnd, thetaAtEnd) };

		// record difference
		double const resid{ gotRefDevAngle - expRefDevAngle };
		resids.emplace_back(resid);

		// check consistency
		std::ostringstream msg;
		constexpr double tolAngle{ .000005 }; // about 1 arc second
		checkRefDevAngles
			( msg
			, gotRefDevAngle
			, expRefDevAngle
			, highSensor
			, highGround
			, tolAngle
			);
		if (! msg.str().empty())
		{
			// append any error case message to test error stream
			oss << msg.str() << '\n';
		}

		// generate residual report for potential later use
		using engabra::g3::io::fixed;
		rptResid
			<< "  Alt[m]: " << fixed(altDev.first, 5u, 0u)
			<< "  Deviation:MoP[uRad]: " << fixed(1.e6*altDev.second, 3u, 1u)
			<< "  Residual(got-exp)[uRad]: " << fixed(1.e6*resid, 3u, 1u)
			<< '\n';
	}

	// summary information
	std::ostringstream rpt;
	rpt << "\n#===Residuals: MoP Validation Results";
	rpt << "\n#===  [MoP: Manual of Photogrammetry (3rd ed., pg 487)]";
	rpt << "\n#===  [45-deg look angle (from Nadir)]";
	rpt << "\n#===  [ground elevation 0.]";
	rpt << '\n' << rptResid.str();
	rpt << "#===\n";
	if (showResiduals)
	{
		std::cout << rpt.str() << std::endl;
	}

	if (! oss.str().empty())
	{
		oss << rptResid.str() << std::endl;
	}

// ExampleEnd

	return oss.str();
}

}

/*! \brief Unit test for Refraction computation
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);
	test1(oss);
	test2(oss);

	return tst::finish(oss);
}


