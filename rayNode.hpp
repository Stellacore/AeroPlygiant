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

#ifndef aply_ray_Node_INCL_
#define aply_ray_Node_INCL_

/*! \file
 *
 * \brief Ray propagation simulation functions.
 *
 */

#include "rayDirChange.hpp"

#include <Engabra>

#include <sstream>
#include <string>


namespace aply
{
namespace ray
{
	using namespace engabra::g3;

	//! Data relevant to an individual node
	struct Node
	{
		Vector const thePrevTan;
		double const thePrevNu;
		Vector const theCurrLoc;
		double const theNextNu;
		Vector const theNextTan;
		DirChange const theDirChange;

		//! Descriptive information about this instance
		inline
		std::string
		infoBrief // Node::
			( std::string const & title = {}
			, std::size_t const & precisionVec = 6u
			, std::size_t const & precisionNu = 6u
			) const
		{
			std::ostringstream oss;
			if (! title.empty())
			{
				oss << title << '\n';
			}
			oss << " tan ";
			oss << io::fixed(thePrevTan, 3u, precisionVec);
			oss << " nu ";
			oss << io::fixed(thePrevNu, 3u, precisionNu);
			oss << " loc ";
			oss << io::fixed(theCurrLoc, 3u, precisionVec);
			oss << " nu ";
			oss << io::fixed(theNextNu, 3u, precisionNu);
			oss << " tan ";
			oss << io::fixed(theNextTan, 3u, precisionVec);
			//oss << " chng ";
			oss << "  ";
			oss << nameFor(theDirChange);
			return oss.str();
		}

		//! Descriptive information about this instance
		inline
		std::string
		infoString // Node::
			( std::string const & title = {}
			) const
		{
			std::ostringstream oss;
			if (! title.empty())
			{
				oss << title << '\n';
			}
			oss << "thePrevTan: " << io::fixed(thePrevTan, 8u, 6u);
			oss << '\n';
			oss << " thePrevNu: " << io::fixed(thePrevNu, 8u, 6u);
			oss << '\n';
			oss << "theCurrLoc: " << io::fixed(theCurrLoc, 8u, 6u);
			oss << '\n';
			oss << " theNextNu: " << io::fixed(theNextNu, 8u, 6u);
			oss << '\n';
			oss << "theNextTan: " << io::fixed(theNextTan, 8u, 6u);
			oss << '\n';
			oss << "theDirChange" << nameFor(theDirChange);
			return oss.str();
		}

		//! Node associated with reversing direction of propagation
		inline
		Node
		reversed // Node::
			() const
		{
			return Node
				{ -theNextTan
				,  theNextNu
				,  theCurrLoc
				,  thePrevNu
				, -thePrevTan
				,  reverseChange(theDirChange)
				};
		}

	}; // Node

} // [ray]
} // [aply]


#endif // aply_ray_Node_INCL_

