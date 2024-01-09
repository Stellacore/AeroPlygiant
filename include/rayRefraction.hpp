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

#ifndef aply_ray_Refraction_INCL_
#define aply_ray_Refraction_INCL_

/*! \file
\brief Declarations for ray::Refraction
*/


#include "envAtmosphere.hpp"

#include <Engabra>

#include <string>
#include <vector>


namespace aply
{
namespace ray
{

/*! \brief Determine angle and displacement due to refraction.

			TODO

\par Example
\dontinclude testgeo/uRefraction.cpp
\skip ExampleStart
\until ExampleEnd

This class computes the displacement that results when a ray of light
(as between the ground and a sensor) travels through the atmosphere.
It uses the COESA1976 atmosphere parameters to obtain the needed
index of refraction data.

It is possible to determine the angular displacement from the line between
the sensor and the center of the earth and also to determine the displacement
along the ground.
*/

class Refraction
{

	//! refraction system
//	friend class RefractionSystem;

public: // methods

	//! default null constructor
	Refraction
		();

	/*! \brief Construct refraction engine to propagate ray.
	 *
	 * Propagation is performed in a local coordinate frame for which:
	 *
	 *	- e3 ('z') axis is direction from Earth center directed
	 *	  vertically upward through sensor station location.
	 *	- Planar coordinate system containing e3 axis and ray path.
	 *	- Earth Radius is distance from center of Earth to ground
	 *	  nadir point.
	 *		- NOTE: computation is not very sensitive to this values
	 *		  so that any reasonable approximation is good enough.
	 *
	 * Computation
	 *	- arg Ray path begins by leaving from radSensor (on e3 axis)
	 *	  above ground nadir location.
	 *
	an instance out of the various parameters needed. radSensor is
	    the radius of the starting position, and startAngle is the angle between
	    the ray and the vertical.
	*/
	explicit
	Refraction
		( double const & lookAngle
		, double const & radSensor
		, double const & radiusEarth
		);

	// destructor -- compiler provided

	//! Check if instance is valid
	bool
	isValid
		() const;

	/*! \brief Determine angle 'Theta_c' deviation.
	 *
	 * Theata is the angle subtended from center of earth between
	 * direction to sensor location and direction to point along ray.
	*/
	double
	thetaAngleAt
		( double const & radius
		) const;

	/*  Determine the displacement from the vertical at a certain
		radius.
	*/
	double
	displacementAt
		( double const & radius
		) const;

	//! Descriptive information about this instance.
	std::string
	infoString
		( std::string const & title=std::string()
		, std::string const & fmt="%20.15g"
		) const;


// TODO private: // data

	env::Atmosphere theAtmosphere{};
	double theRadiusEarth{ engabra::g3::null<double>() };
	double theRefractiveInvariant{ engabra::g3::null<double>() };
	mutable std::pair<double, std::vector<double> > theInitValues{};
};

} // ray
} // aply

#endif // aply_ray_Refraction_INCL_
