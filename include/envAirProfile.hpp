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

#ifndef aply_env_AirProfile_INCL_
#define aply_env_AirProfile_INCL_

/*! \file
\brief Declarations for aply::env::AirProfile
*/


#include "envAirInfo.hpp"


namespace aply
{
namespace env
{

/*! \brief Wrapper to interpolate AirInfo data from (ordered) collections.
*/

struct AirProfile
{
	//! Collection of AirInfo properties ordered by height above ground.
	std::map<Height, AirInfo> theAirInfoMap{};

	//! AirInfo values interpolated at given high above ground (elevation).
	AirInfo
	airInfoAtHeight
		( double const & height
		) const;

	//! Index of refraction corresponding to this height.
	double
	indexOfRefraction
		( double const & height
		) const;

	//! True if theAirInfoMap has at least two entries (needed to interpolate).
	bool
	isValid
		() const;
};

} // [env]
} // [aply]

#endif //  aply_env_AirProfile_INCL_

