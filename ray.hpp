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

#ifndef aply_ray_INCL_
#define aply_ray_INCL_

/*! \file
 *
 * \brief Ray propagation simulation functions.
 *
 */


#include "rayDirChange.hpp"
#include "rayNode.hpp"
#include "rayPath.hpp"
#include "rayPathView.hpp"
#include "rayPropagator.hpp"
#include "rayStart.hpp"

#include <iostream>


namespace aply
{
/*! \brief Functions and classes for simulation of ray propagation
 */
namespace ray
{

} // [ray]
} // [aply]


namespace
{

	//! Overload output for ray::Start
	inline
	std::ostream &
	operator<<
		( std::ostream & ostrm
		, aply::ray::Start const & start
		)
	{
		ostrm << start.infoString();
		return ostrm;
	}

	//! Overload output for ray::Node
	inline
	std::ostream &
	operator<<
		( std::ostream & ostrm
		, aply::ray::Node const & node
		)
	{
		ostrm << node.infoString();
		return ostrm;
	}

} // [anon]


#endif // aply_ray_INCL_

