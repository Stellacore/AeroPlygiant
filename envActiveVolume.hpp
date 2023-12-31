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

#ifndef aply_env_ActiveVolume_INCL_
#define aply_env_ActiveVolume_INCL_

/*! \file
 *
 * \brief Environment configuration parameters (related to Refraction)
 *
 */


#include <Engabra>

#include <memory>
#include <string>



namespace aply
{
namespace env
{
	using namespace engabra::g3;

	/*! \brief Specify volume of space through which rays should be propagated.
	 */
	struct ActiveVolume
	{
		std::string theName{};

		//! Construct a named instance
		explicit
		ActiveVolume // ActiveVolume::
			( std::string const & name = "ActiveVolume"
			)
			: theName{ name }
		{ }

		//! Overload to define shape of volume (true: inside, false: outside)
		inline
		virtual
		bool
		contains // ActiveVolume::
			( Vector const & rVec
			) const
		{
			return true;
		}

	}; // ActiveVolume

	//! An active volume w/o limits
	static std::shared_ptr<ActiveVolume> const sPtAllSpace
		{ std::make_shared<ActiveVolume>("sAllSpace") };

	//! A rectangular ActiveVolume determined by two corner points.
	struct ActiveBox : public ActiveVolume
	{
		//! True if (minIncluded <= value < maxExcluded)
		inline
		static
		bool
		inInterval
			( double const & minIncluded
			, double const & value
			, double const & maxExcluded
			)
		{
			return ( (! (value < minIncluded)) && (value < maxExcluded) );
		}

		Vector const theMinCorner{ null<Vector>() };
		Vector const theMaxCorner{ null<Vector>() };

		inline
		explicit
		ActiveBox
			( Vector const & minCorner
			, Vector const & maxCorner
			)
			: theMinCorner{ minCorner }
			, theMaxCorner{ maxCorner }
		{ }

		//! true if theMinCorner[ndx] <= rVec[ndx] < theMaxCorner[ndx] all ndx.
		inline
		virtual
		bool
		contains
			( Vector const & rVec
			) const
		{
			return
				(  inInterval(theMinCorner[0], rVec[0], theMaxCorner[0])
				&& inInterval(theMinCorner[1], rVec[1], theMaxCorner[1])
				&& inInterval(theMinCorner[2], rVec[2], theMaxCorner[2])
				);
		}


	}; // ActiveBox

} // [env]
} // [aply]

#endif // aply_env_ActiveVolume_INCL_

