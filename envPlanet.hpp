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

#ifndef aply_env_Planet_INCL_
#define aply_env_Planet_INCL_

/*! \file
 *
 * \brief Environment configuration parameters (related to Refraction)
 *
 */


#include <Engabra>

#include <sstream>
#include <string>


namespace aply
{
namespace env
{
	using namespace engabra::g3;

	/*! \brief Data parameters for spherical planet's atmosphere.
	 *
	 * Container for data describing useful for creating a simple
	 * spherically symmetric index of refraction field (e.g. for
	 * describing Earth atmospher models)
	 *
	 * The 'nu' naming refer to Index of Refraction (IoR) values.
	 */
	struct Planet
	{

		//! Index of refraction ('nu') at ground level.
		double const theNuGround{ null<double>() };
		//! Index of refraction ('nu') at edge of space
		double const theNuSpace{ null<double>() };

		//! Radius from planet center with IoR == theNuGround
		double const theRadGround{ null<double>() };
		//! Radius from planet center with IoR == theNuSpace
		double const theRadSpace{ null<double>() };

		//! Descriptive information about this instance
		inline
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
			oss << " theNuGround: " << io::fixed(theNuGround);
			oss << '\n';
			oss << "  theNuSpace: " << io::fixed(theNuSpace);
			oss << '\n';
			oss << "theRadGround: " << io::fixed(theRadGround);
			oss << '\n';
			oss << " theRadSpace: " << io::fixed(theRadSpace);
			return oss.str();
		}

	}; // Planet

	//! Parameters for Earth atmosphere
	static Planet const sEarth
		{
		  1.000273 // nuGround // (for earth, actual is 1.000273 at STP)
		, 1.000    // nuSpace
		, 6370.e3  // radiusGround
		, 6470.e3  // radiusSpace
		};

} // [env]
} // [aply]

#endif // aply_env_Planet_INCL_

