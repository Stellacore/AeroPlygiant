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
 * \brief Unit test for class env::Air
 *
 */


#include "envAir.hpp"

#include "tst.hpp"

#include <sstream>


/* Example data from: http://weather.uwyo.edu/upperair/sounding.html

>>>
   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
    hPa     m      C      C      %    g/kg    deg   knot     K      K      K 
-----------------------------------------------------------------------------
 1000.0    136                                                               
  925.0    747                                                               
  913.0    849   -8.5  -10.3     87   1.92    325     11  271.6  277.1  271.9
  908.0    892   -9.7  -13.6     73   1.48    317     14  270.8  275.1  271.1
  907.0    900   -9.8  -13.6     74   1.48    315     15  270.8  275.1  271.1
<<<

*/

namespace
{


	//! Check basic env::AirInfo methods
	void
	test0
		( std::ostringstream & oss
		)
	{
		// [DoxyExample00]
		// [DoxyExample00]

		using namespace aply::env;
		// Example fields from UWyo format (some spaces elided)
		constexpr char sUWyoLine[] = 
			"  908.0  892 -9.7 -13.6  73 1.48   317   14  270.8  275.1  271.1";
		constexpr double expHigh{ 892. };
		constexpr double expTemp{ -9.7 + 273.15 };
		constexpr double expPres{ 908.0 * 100. };
		constexpr double expRelH{ 73. * .01 };
		AirInfo const info{ AirInfo::fromUWyoRecord(sUWyoLine) };
		double const & gotHigh = info.theHigh;
		double const & gotTemp = info.theTemp;
		double const & gotPres = info.thePres;
		double const & gotRelH = info.theRelH;

		using engabra::g3::nearlyEquals;
		if (! nearlyEquals(gotHigh, expHigh))
		{
			oss << "Failure of High test\n";
			oss << "expHigh: " << expHigh << '\n';
			oss << "gotHigh: " << gotHigh << '\n';
		}
		if (! nearlyEquals(gotTemp, expTemp))
		{
			oss << "Failure of Temp test\n";
			oss << "expTemp: " << expTemp << '\n';
			oss << "gotTemp: " << gotTemp << '\n';
		}
		if (! nearlyEquals(gotPres, expPres))
		{
			oss << "Failure of Pres test\n";
			oss << "expPres: " << expPres << '\n';
			oss << "gotPres: " << gotPres << '\n';
		}
		if (! nearlyEquals(gotRelH, expRelH))
		{
			oss << "Failure of RelH test\n";
			oss << "expRelH: " << expRelH << '\n';
			oss << "gotRelH: " << gotRelH << '\n';
		}

std::cout << "uwyoLine: " << info << std::endl;

	}
}


/*! \brief Unit test for env::Air
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);

	return tst::finish(oss);
}

