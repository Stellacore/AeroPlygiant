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
\brief Definitions for aply::env::AtmosphereParameters
*/


#include "envAtmosphereParameters.hpp"

#include <Engabra>

#include <sstream>


namespace aply
{
namespace env
{

bool
AtmosphereParameters :: isValid() const
{
	return engabra::g3::isValid(theHigh)
		&& engabra::g3::isValid(theTemp)
		&& engabra::g3::isValid(thePressure)
		&& engabra::g3::isValid(theIoR);
}

std::string
AtmosphereParameters :: infoBrief
	( std::string const & title
	) const
{
	std::ostringstream oss;
	if (! title.empty())
	{
		oss << std::setw(12u) << title << ' ';
	}
	using engabra::g3::io::fixed;
	oss
		<< " H[m],T[K],P[mBar],IoR[-]: "
		<< " " << fixed(theHigh, 4u, 3u)
		<< " " << fixed(theTemp, 4u, 2u)
		<< " " << fixed(thePressure, 4u, 1u)
		<< " " << fixed(theIoR, 1u, 9u)
		;
	return oss.str();
}

std::string
AtmosphereParameters :: infoString
	( std::string const & title
	) const
{
	std::ostringstream oss;

	if (! title.empty())
	{
		oss << title << std::endl;
	}

	using engabra::g3::io::fixed;
	oss << "Height:   " << fixed(theHigh) << std::endl;
	oss << "Temp:     " << fixed(theTemp) << std::endl;
	oss << "Pressure: " << fixed(thePressure) << std::endl;
	oss << "RefIndex: " << fixed(theIoR);

	return oss.str();
}

} // [env]
} // [aply]

