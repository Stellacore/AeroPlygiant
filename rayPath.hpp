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

#include "rayStart.hpp"

//#include "env.hpp"
#include "rayNode.hpp"

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
	 * Monitors path propagation progress in order to determine
	 * when to stop propagation computations. Intermediate results
	 * are archived (into theNodes) along the way.
	 *
	 * After consuming propagation data, theNodes collection, provides
	 * a representation of the overall curved (polygonal) propagation
	 * path.
	 *
	 * Provides methods that are compatible with those of std::vector.
	 * The capacity() is originally set based on nominal distance from
	 * the theStart location. After which capacity() is reported larger
	 * than the size() while the ray is closing on theStopLoc. When the
	 * ray begins opening distance to theStopLoc, the capacity() is
	 * set to zero (which forces the Propagator to stop.
	 *
	 * TODO - need a smarter way to control propagation length (or at
	 *        least a better way to implement the idea).
	 *
	 * TODO - Factor shape/info reporting (wrappers to use collection of Nodes)
	 *
	 */
	struct Path
	{
		//! Starting boundary condition (direction and location) for the ray
		Start const theStart{};
		//! Point used for estimating path size used for Propagation
		Vector const theStopLoc{ null<Vector>() };
		//! Increment specifying how often to archive path data in theNodes.
		double const theSaveDelta{ null<double>() };
		//! Archived path information (approximately every theSaveDelta units)
		std::vector<ray::Node> theNodes{};

	private:

		//! Tracking value (how close to theStopLoc on previous emplace_back()
		double thePrevNearDist{ null<double>() };
		//! Tracking value (how close to theStopLoc currently)
		double theCurrNearDist{ null<double>() };

	public:

		//! Construct storage based on nominal distance between points
		inline
		explicit
		Path // Path::
			( Start const & startWith
				//!< Initial direction and start point for propagation
			, Vector const & stopNearTo
				//!< Stop tracing with path stops getting closer to this
			, double const & saveStepSize
				//!< Save node if path exceeds this distance from previous save
			)
			: theStart{ startWith }
			, theStopLoc{ stopNearTo }
			, theSaveDelta{ saveStepSize }
			, thePrevNearDist{ null<double>() }
			, theCurrNearDist{ 1.e10 }
		{
			// estimate distance (as if straight line)
			if (engabra::g3::isValid(theStopLoc))
			{
				theCurrNearDist = magnitude(theStopLoc - theStart.thePntLoc);
				double const nomDist
					{ magnitude(theStopLoc - theStart.thePntLoc) };
				constexpr double padFactor{ 9./8. }; // about 12% extra
				double const dubSize{ padFactor * nomDist / theSaveDelta };
				std::size_t const nomSize{ static_cast<std::size_t>(dubSize) };
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

		//! Indicate how much this instance *can* store.
		inline
		std::size_t
		capacity // Path::
			() const
		{

			std::size_t cap{ 0u }; // default is to stop
			if (keepGoing())
			{
				cap = theNodes.capacity();
			}
			return cap;
		}

		//! Process a node - determine if should be archived or not
		inline
		void
		emplace_back // Path::
			( ray::Node const & node
			)
		{
			updateNearDists(node);
			double distFromSave{ 0. };
			if (! theNodes.empty())
			{
				Vector const pathDelta
					{ node.theCurrLoc - theNodes.back().theCurrLoc };
				distFromSave = magnitude(pathDelta);
			}

			bool const isFirstStep{ theNodes.empty() };
			bool const pastStepSize{ ! (distFromSave < theSaveDelta) };
			bool const aboutToStop{ ! keepGoing() };

			/*
			std::cout
				<< "  thePrevNearDist: " << io::fixed(thePrevNearDist) << '\n'
				<< "  theCurrNearDist: " << io::fixed(theCurrNearDist) << '\n'
				<< "  theCurrLoc: " << io::fixed(node.theCurrLoc) << '\n'
				<< "  theStopLoc: " << io::fixed(theStopLoc) << '\n'
				<< "     distFromSave: " << io::fixed(distFromSave) << '\n'
				<< "  T/F:isFirstStep: " << isFirstStep << '\n'
				<< " T/F:pastStepSize: " << pastStepSize << '\n'
				<< "  T/F:aboutToStop: " << aboutToStop << '\n'
				<< "             size: " << size() << '\n'
				<< "         capacity: " << capacity() << '\n'
				;
			*/

			if (isFirstStep || pastStepSize || aboutToStop)
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
			oss << "theStopLoc: " << theStopLoc;
			oss << '\n';
			oss << "theSaveDelta: " << theSaveDelta;
			oss << '\n';
			oss << "theNodes.size(): " << theNodes.size();
		//	oss << '\n';
		//	oss << "thePrevNearDist: " << io::fixed(thePrevNearDist);
		//	oss << '\n';
		//	oss << "theCurrNearDist: " << io::fixed(theCurrNearDist);

			return oss.str();
		}

	private:

		//! True when the nearest distance to stop point starts increasing
		inline
		bool
		keepGoing // Path::
			() const
		{
			bool keepGo{ true };
			if (isValid(theCurrNearDist) && isValid(thePrevNearDist))
			{
				keepGo = ! (thePrevNearDist < theCurrNearDist);
			}
			return keepGo;
		}

		//! Archive this node in storage
		inline
		void
		addNode // Path::
			( ray::Node const & node
			)
		{
			theNodes.emplace_back(node);
		}

		//! Distance (cord length - NOT along ray) to stop request point
		inline
		double
		distanceFromStop // Path::
			( ray::Node const & node
			)
		{
			return magnitude(theStopLoc - node.theCurrLoc);
		}

		//! Keep track of Previous and Current nearest distances
		inline
		void
		updateNearDists // Path::
			( ray::Node const & currNode
			)
		{
			// compute the current distance to stopping point
			double const currDist{ distanceFromStop(currNode) };
			// ping pong the buffers
			thePrevNearDist = theCurrNearDist;
			theCurrNearDist = currDist;
		}
	
	}; // Path

} // [ray]
} // [aply]


#endif // aply_ray_Path_INCL_

