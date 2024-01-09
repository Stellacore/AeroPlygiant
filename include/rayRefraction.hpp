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

This class computes the displacement that results when a ray of light
travels through the atmosphere.  (as between a sensor and a location on
the ground).

It uses the COESA1976 atmosphere parameters to obtain the needed
index of refraction data.

The ray path is expressed in terms of polar coordinates relative to
origin at the \b center \b of \b Earth. I.e. a point on the ray path has
location (r,theta) where "r" is a value on the order of 6.37e6 [m] and
theta=0 is start point of ray.

The ray initial conditions (start location and direction) are provided
to the constructor. The construction is light weight (e.g. nothing is
compuated initially).

At any point(s) in the future, the constructed instance may be queried
to obtain an end point on the ray - by using thetaAngleAt() method.

The thetaAngleAt() method, performs the full numeric integration
computations (e.g. is the potentially "expensive" operation). Note
that this is the polar angle (from center of Earth) at which the
ray passes through a distance (from center of Earth) equal to the
radiusEnd argument to the thetaAngleAt() function.

The deflection angle observed \b from the \b sensor position needs
to be computed using the initial conditions (location and direction)
and the obtained end point. E.g.

\snippet test/test_Refraction.cpp DoxyExample01

The refraction model is that presented by Gyer 1996:
\verbatim
@article{gyer1996:AtmRefraction,
	title = {Methods for Computing Photogrammetric Refraction Corrections for Vertical and Oblique Photographs},
	author = {Maurice S. Gyer},
	journal = {Photogrammetric Engineering \& Remote Sensing},
	year = {1996},
	month = {March},
	pages = {301-310},
	url = {https://www.asprs.org/wp-content/uploads/pers/1996journal/mar/1996_mar_301-310.pdf},
	urldate = {2023-12-08},
}
\endverbatim

*/

class Refraction
{

private: // data

	//! Cached construction value.
	double const theStartLookAngle{ engabra::g3::null<double>() };

	//! Cached construction value.
	double const theStartRadius{ engabra::g3::null<double>() };

	/*! \brief Initial value of Theta_c.
	 *
	 * Here this is always zero since ray tracing is performed in a
	 * local coordinate system (e.g. the polar axis is constructed to 
	 * pass through the sensor station by definition).
	 */
	static constexpr double theTheta0{ 0. };

	/*! \brief Atmospheric model providing IoR as function of \b elevation.
	 *
	 * \note The atmospheric model is queried for an IoR value at a
	 * \b elevation value (e.g. at a current height above the radiusEarth
	 * value provided to Refraction() ctor. The reason for doing this
	 * is that it keeps the atmospheric model relatively DE-coupled
	 * from any specific figure of Earth models.
	 */
	env::Atmosphere theAtmosphere{};

	//! \brief Defines the "zero-elevation" location relative to ECEF origin.
	double theRadiusEarth{ engabra::g3::null<double>() };

	/*! \brief The Snell's constant (IoR * sin(angle)) value.
	 *
	 * This is the "k" value in Gyer's paper, Eqn[1].
	 */
	double theRefractiveInvariant{ engabra::g3::null<double>() };

	/*! \brief Initial conditions - polar coordinate of ray starting location.
	 *
	 * Initial value structure includes:
	 *	- theInitRadTheta.first : radius (should match method input argument)
	 *	- theInitRadTheta.second: has size of 1u and contains
	 *		- [0] : Theta_c angle (ray path polar angle from center of Earth)
	 */
	std::pair<double, std::vector<double> > theInitRadTheta{};

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
	 * Computation includes:
	 *	- Ray path begins by leaving sensor from radiusSensor distance
	 *	  from Earth center (on e3 axis). Ref ctor.
	 *	- Ray path leaves sensor in lookAngle (from Nadir) direction
	 *	  (e.g. 0 angle is ray heading straight down). Ref ctor.
	 *	- Ray propagates until distance radiusEnd from Earth center.
	 *	  Ref thetaAngleAt() and displacementAt().
	*/
	explicit
	Refraction
		( double const & lookAngle
		, double const & radiusSensor
		, double const & radiusEarth
		);

	// destructor -- compiler provided

	//! Check if instance is valid
	bool
	isValid
		() const;

	/*! \brief Determine angle 'Theta_c' deviation.
	 *
	 * Theata is the angle subtended from center of Earth between
	 * direction to sensor location (i.e. the positive e3 axis) and
	 * the direction (from Earth center) to (end) point of ray a
	 * distance radiusEnd from Earth center.
	 *
	 * Ref Gyer 1996 Fig 2.
	 *
	 * \note This function performs numerical integration computations
	 *       and therefore can take a non-trivial amount of time.
	 */
	double
	thetaAngleAt
		( double const & radiusEnd
		) const;

	/*! \brief Angular deviation of ray end as observed from start point.
	 *
	 * The ray leaves the start point (ref #theInitRadTheta) at a look
	 * angle (relative to Nadir direction). Call this the "observed look
	 * angle".
	 *
	 * The ray follows a curved path which terminates at the "EndPoint"
	 * (specified by combination of radiusEnd and thetaAngleAt(radiusEnd)).
	 *
	 * From the Start point a geometrically straight line toward the
	 * EndPoint defines the "ideal look angle".
	 *
	 * This function returns "ideal look angle" minus "observed look angle".
	 */
	double
	angularDeviationFromStart
		( double const radiusEnd
			//!< Distance from \b center of Earth at which ray terminates.
		, double const thetaEnd
			//!< Angle subtended by ray path \b from \b Earth \b center.
		) const;

	/*! \brief Convienience: angularDeviationFromStart(thetaAngleAt(radiusEnd))
	 *
	 * \note: This performs numeric integration and therefore may take
	 *        non trivial amount of time.
	 */
	inline
	double
	angularDeviationFromStart
		( double const radiusEnd
			//!< Distance from \b center of Earth at which ray terminates.
		) const
	{
		return angularDeviationFromStart(radiusEnd, thetaAngleAt(radiusEnd));
	}

	//! Descriptive information about this instance.
	std::string
	infoString
		( std::string const & title=std::string()
		, std::string const & fmt="%20.15g"
		) const;

};

} // ray
} // aply

#endif // aply_ray_Refraction_INCL_
