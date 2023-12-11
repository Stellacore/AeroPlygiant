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

#ifndef aply_ray_Path_INCL_
#define aply_ray_Path_INCL_

/*! \file
 *
 * \brief Ray propagation simulation functions.
 *
 */


#include "rayNode.hpp"
#include "rayStart.hpp"

#include <Engabra>

#include <sstream>
#include <string>
#include <vector>


namespace aply
{
namespace ray
{
	using namespace engabra::g3;

	/*! \brief Consumer of dynamically generated path information.
	 *
	 * Wraps a collection (std::vector) of Node instances. The
	 * considerNode() method monitors path propagation length since
	 * the previous Node was added to the collection. If the pathlength
	 * has increased by more then #theSaveDelta amount, then the
	 * considered node is added to #theNodes collection.
	 * 
	 * Provides methods that are compatible with those of std::vector
	 * so that an instance can be used in generic programming contexts.
	 */
	struct Path
	{
		//! Starting boundary condition (direction and location) for the ray
		Start const theStart{};
		//! Increment specifying how often to archive path data in theNodes.
		double const theSaveDelta{ null<double>() };
		//! Archived path information (approximately every theSaveDelta units)
		std::vector<ray::Node> theNodes{};

		//! Estimate collection size needed to span between beg/end locations.
		inline
		static
		std::size_t
		sizeBetween
			( Vector const & begLoc
			, Vector const & endLoc
			, double const & deltaDist
			, double const & padFactor = 9./8.
			)
		{
			// estimate distance (as if straight line)
			double const nomDist{ magnitude(endLoc - begLoc) };
			// make a bit larger to allow for path curvature/changes
			double const dubSize{ padFactor * nomDist / deltaDist };
			std::size_t const nomSize{ static_cast<std::size_t>(dubSize) };
			return nomSize;
		}


		//! Construct storage based on nominal distance between points
		inline
		explicit
		Path // Path::
			( Start const & startWith
				//!< Initial direction and start point for propagation
			, double const & saveStepSize
				//!< Save node if path exceeds this distance from previous save
			, Vector const & approxEndLoc = null<Vector>()
				//!< Used to estimate/allocate storage space
			)
			: theStart{ startWith }
			, theSaveDelta{ saveStepSize }
		{
			// estimate distance (as if straight line)
			if (engabra::g3::isValid(approxEndLoc))
			{
				Vector const & begLoc = theStart.thePntLoc;
				std::size_t const nomSize
					{ sizeBetween(begLoc, approxEndLoc, theSaveDelta) };
				theNodes.reserve(nomSize);
			}
		}

		//! Indicate how much this instance currently *has* stored.
		inline
		std::size_t
		size // Path::
			() const
		{
			return theNodes.size();
		}

		/*! \brief Set maximum capacity.
		 * 
		 * This value controls termination for cases where the
		 * media has valid IoR value for long (or infinite) distances.
		 */
		inline
		void
		reserve
			( std::size_t const & maxNodeSize
			)
		{
			theNodes.reserve(maxNodeSize);
		}

		//! Reserve enough space for this (arc-length) at #theSaveDelta.
		inline
		void
		reserveForDistance
			( double const & dist
			)
		{
			if (theSaveDelta < dist)
			{
				double const dNum{ dist / theSaveDelta };
				std::size_t const numElem{ static_cast<std::size_t>(dNum) + 1u};
				theNodes.reserve(numElem);
			}
		}

		//! Indicate how much this instance *can* store.
		inline
		std::size_t
		capacity // Path::
			() const
		{
			return theNodes.capacity();
		}

		//! Process a node - determine if should be archived or not
		inline
		void
		emplace_back // Path::
			( ray::Node const & node
			)
		{
			considerNode(node);
		}

		//! Process a node - determine if should be archived or not
		inline
		void
		considerNode // Path::
			( ray::Node const & node
			)
		{
			double distFromSave{ 0. };
			if (! theNodes.empty())
			{
				Vector const pathDelta
					{ node.theCurrLoc - theNodes.back().theCurrLoc };
				distFromSave = magnitude(pathDelta);
			}

			bool const isFirstStep{ theNodes.empty() };
			bool const pastStepSize{ ! (distFromSave < theSaveDelta) };

			if (isFirstStep || pastStepSize)
			{
				addNode(node);
			}
		}

		//! Descriptive information about this instance
		inline
		std::string
		infoString // Path::
			( std::string const & title = {}
			) const
		{
			std::ostringstream oss;
			if (! title.empty())
			{
				oss << title << '\n';
			}
			oss << "theStart: " << theStart.infoString();
			oss << '\n';
			oss << "theSaveDelta: " << theSaveDelta;
			oss << '\n';
			oss << "theNodes.size(): " << theNodes.size()
				<< "  of(capacity)  " << theNodes.capacity();

			return oss.str();
		}

	private:

		//! Archive this node in storage
		inline
		void
		addNode // Path::
			( ray::Node const & node
			)
		{
			theNodes.emplace_back(node);
		}
	
	}; // Path

} // [ray]
} // [aply]


#endif // aply_ray_Path_INCL_

