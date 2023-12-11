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
	 * has increased by more then #theSaveDist amount, then the
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
		double const theSaveDist{ null<double>() };
		//! Archived path information (approximately every theSaveDist units)
		std::vector<ray::Node> theNodes{};

	private:

		//! Track (approximate) residual arc-length since last archived node
		double theResidArcDist{ null<double>() };
		//! The location of the last considered (but generally not saved) node
		Vector theLastSeenLoc{ null<Vector>() };

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

	public:

		//! Construct storage based on nominal distance between points
		inline
		explicit
		Path // Path::
			( Start const & startWith
				//!< Initial direction and start point for propagation
			, double const & saveStepDist
				//!< Save node if path exceeds this distance from previous save
			, Vector const & approxEndLoc = null<Vector>()
				//!< Used to estimate/allocate storage space
			)
			: theStart{ startWith }
			, theSaveDist{ saveStepDist }
			, theNodes{ }
			, theResidArcDist{ 0. }
			, theLastSeenLoc{ null<Vector>() }
		{
			// estimate distance (as if straight line)
			if (engabra::g3::isValid(approxEndLoc))
			{
				Vector const & begLoc = theStart.thePntLoc;
				std::size_t const nomSize
					{ sizeBetween(begLoc, approxEndLoc, theSaveDist) };
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

		//! Reserve enough space for this (arc-length) at #theSaveDist.
		inline
		void
		reserveForDistance
			( double const & dist
			)
		{
			if (theSaveDist < dist)
			{
				double const dNum{ dist / theSaveDist };
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
			bool saveThisNode{ false };

			if (theNodes.empty())
			{
				saveThisNode = true;
			}
			else
			{
				// check distance from previously considered node
				Vector const delta{ node.theCurrLoc - theLastSeenLoc };
				double const deltaMag{ magnitude(delta) };
				// increment residual arc length by this much
				theResidArcDist += deltaMag;
			}

			// check if save distance is exceeded
			if (! (theResidArcDist < theSaveDist))
			{
				saveThisNode = true;
			}

			if (saveThisNode)
			{
				// always add the first node
				theNodes.emplace_back(node);
				// set residual arc distance
				theResidArcDist = 0.;
			}

			// remember the last considered node (whether added or not)
			theLastSeenLoc = node.theCurrLoc;
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
			oss << "theSaveDist: " << theSaveDist;
			oss << '\n';
			oss << "theNodes.size(): " << theNodes.size()
				<< "  of(capacity)  " << theNodes.capacity();

			return oss.str();
		}

	}; // Path

} // [ray]
} // [aply]


#endif // aply_ray_Path_INCL_

