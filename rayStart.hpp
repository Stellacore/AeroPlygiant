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

#ifndef aply_ray_Start_INCL_
#define aply_ray_Start_INCL_

/*! \file
 *
 * \brief Ray propagation simulation functions.
 *
 */


#include <Engabra>

#include <sstream>
#include <string>


namespace aply
{
namespace ray
{
	using namespace engabra::g3;

	/*! \brief Data representing initial boundary values for a ray(curve).
	 *
	 */
	struct Start
	{
		Vector const theTanDir{}; //!< Incident tangent direction (unitary)
		Vector const thePntLoc{}; //!< Point of incidence for tangent dir

		//! Create an instance ensuring tangent dir is unitary.
		inline
		static
		Start
		from // Start::
			( Vector const & anyTan
			, Vector const & loc
			)
		{
			return { direction(anyTan), loc };
		}

		//! Descriptive information about this instance
		inline
		std::string
		infoString // Start::
			( std::string const & title = {}
			) const
		{
			std::ostringstream oss;
			if (! title.empty())
			{
				oss << title << " ";
			}
			oss
				<< "dir: " << theTanDir
				<< ' '
				<< "loc: " << thePntLoc
				;
			return oss.str();
		}

	}; // Start

} // [ray]
} // [aply]


#endif // aply_ray_Start_INCL_

