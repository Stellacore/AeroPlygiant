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

#ifndef aply_env_Atmosphere_INCL_
#define aply_env_Atmosphere_INCL_

/*! \file
\brief Declarations for aply::env::Atmosphere
*/

#include "envAtmosphereParameters.hpp"

#include <string>


namespace aply
{
namespace env
{

/*! \brief Provide estimates of atmospheric data from given atmospheric data.

This algorithm uses a set of data points and linearly interpolates them.

\par Example
\dontinclude test/uAtmosphere.cpp
\skip ExampleStart
\until ExampleEnd
*/

class Atmosphere
{

public: // data

	std::map<double, AtmosphereParameters> theParms;

public: // static

	//! Use COESA1976 model
	static
	Atmosphere
	COESA1976
		();

public: // methods

	//! default null constructor
	Atmosphere
		();

	//! Interpolate all values for given high above sea level
	AtmosphereParameters
	parametersForHeight
		( double const & high
		) const;

// TODO
	inline
	double
	indexOfRefraction
		( double const & pointElevation
		) const
	{
		return -1.; // TODO
	}

	//! Check if instance is valid
	bool
	isValid
		() const;

	//! Descriptive information about this instance.
	std::string
	infoString
		( std::string const & title=std::string()
		, std::string const & fmt = "%20.15g"
		) const;

};

} // [env]
} // [aply]

#endif //  aply_env_Atmosphere_INCL_

