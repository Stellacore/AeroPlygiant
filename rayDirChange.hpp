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

#ifndef aply_ray_DirChange_INCL_
#define aply_ray_DirChange_INCL_

/*! \file
 *
 * \brief Ray propagation simulation functions.
 *
 */

#include <string>


namespace aply
{
namespace ray
{

	//! Characterization of ray path tangent interacting at step boundary.
	enum DirChange
	{
		  Null      //!< Unset or unknown
		, Unaltered //!< Tangent dir unchanged (no gradient)
		, Converged //!< Tangent dir refracted toward gradient (into dense)
		, Diverged  //!< Tangent dir refracted away from gradient (into sparse)
		, Reflected //!< Tangent dir reflected from boundary (total internal)
		, Stopped   //!< Out of simulation domain
		, Started   //!< Begin of ray

	}; // DirChange

	//! Enum value associated with a ray propagating in opposite direction
	inline
	DirChange
	reverseChange
		( DirChange const & fwdChange
		)
	{
		DirChange revChange{ fwdChange }; // Null, Unaltered, Reflected cases
		if (Converged == fwdChange)
		{
			revChange = Diverged;
		}
		else
		if (Diverged == fwdChange)
		{
			revChange = Converged;
		}
		if (Stopped == fwdChange)
		{
			revChange = Started;
		}
		if (Started == fwdChange)
		{
			revChange = Stopped;
		}
		return revChange;
	}

	//! String to associate with each DirChange enum value.
	inline
	std::string
	nameFor
		( DirChange const & change
		)
	{
		std::string name("Unknown");
		switch (change)
		{
			case Unaltered: name = "Unaltered"; break;
			case Converged: name = "Converged"; break;
			case Diverged:  name = "Diverged";  break;
			case Reflected: name = "Reflected"; break;
			default: name = "Null"; break;
		}
		return name;
	}

} // [ray]
} // [aply]


#endif // aply_ray_DirChange_INCL_

