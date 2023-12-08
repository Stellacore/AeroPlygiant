//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//

#ifndef Refraction_save_INCL_
#define Refraction_save_INCL_

/*! \file
 *
 * \brief Environment configuration parameters (related to Refraction)
 *
 */


#include "ray.hpp"


/*! \brief Classes and functions for recording and saving propagated path data.
 */
namespace save
{

	using namespace engabra::g3;

	/*! \brief Store path information suitable for later visualization.
	 *
	 * Basic methods are compatible with those of std::vector.
	 *
	 */
	struct Path
	{
		Vector const theBegTan{ null<Vector>() };
		Vector const theBegLoc{ null<Vector>() };
		Vector const theStopLoc{ null<Vector>() };
		double const theSaveDelta{ null<double>() };
		std::vector<ray::Node> theNodes{};
		double thePrevNearDist{ null<double>() };
		double theCurrNearDist{ null<double>() };

		//! Construct storage based on nominal distance between points
		inline
		explicit
		Path
			( Vector const & beginTangent
				//!< Initial incident direction of propagation
			, Vector const & beginLocation
				//!< Point at which incident ray starts path
			, Vector const & stopNearTo
				//!< Stop tracing with path stops getting closer to this
			, double const & saveStepSize
				//!< Save node if path exceeds this distance from previous save
			)
			: theBegTan{ direction(beginTangent) }
			, theBegLoc{ beginLocation }
			, theStopLoc{ stopNearTo }
			, theSaveDelta{ saveStepSize }
			, thePrevNearDist{ null<double>() }
			, theCurrNearDist{ magnitude(theStopLoc - theBegLoc) }
		{
			// estimate distance (as if straight line)
			double const nomDist{ magnitude(stopNearTo - beginLocation) };
			constexpr double padFactor{ 9./8. }; // about 12% extra
			double const dubSize{ padFactor * nomDist / theSaveDelta };
			std::size_t const nomSize{ static_cast<std::size_t>(dubSize) };
			theNodes.reserve(nomSize);
		}

		//! Indicate how much this instance currently *has* stored.
		inline
		std::size_t
		size
			() const
		{
			return theNodes.size();
		}

		//! True when the nearest distance to stop point starts increasing
		inline
		bool
		keepGoing
			() const
		{
			bool keepGo{ true };
			if (isValid(theCurrNearDist) && isValid(thePrevNearDist))
			{
				keepGo = ! (thePrevNearDist < theCurrNearDist);
			}
			return keepGo;
		}

		//! Indicate how much this instance *can* store.
		inline
		std::size_t
		capacity
			() const
		{
			std::size_t cap{ 0u }; // default is to stop
			if (keepGoing())
			{
				cap = theNodes.capacity();
			}
			return cap;
		}

		//! Archive this node in storage
		inline
		void
		addNode
			( ray::Node const & node
			)
		{
			theNodes.emplace_back(node);
		}

		double
		distanceFromStop
			( ray::Node const & node
			)
		{
			return magnitude(theStopLoc - node.theCurrLoc);
		}

		void
		updateNearDists
			( ray::Node const & currNode
			)
		{
			// compute the current distance to stopping point
			double const currDist{ distanceFromStop(currNode) };
			// ping pong the buffers
			thePrevNearDist = theCurrNearDist;
			theCurrNearDist = currDist;
		}

		//! Process a node - determine if should be archived or not
		inline
		void
		emplace_back
			( ray::Node const & node
			)
		{
//std::cout << node.infoBrief() << '\n';
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
//std::cout << "**** SAVING ****\n";
//std::cout << "\n\n **** SAVING ****\n\n";
				addNode(node);
			}
		}

	}; // Path

} // [save]

#if 0
	/*! \brief Interface specification for consuming structure
	 *
	 * Basic methods are compatible with those of std::vector.
	 *
	 */
	struct Store
	{
		//! Indicate how much this instance currently *has* stored.
		virtual
		std::size_t
		size
			() const = 0;

		//! Indicate how much this instance *can* store.
		virtual
		std::size_t
		capacity
			() const = 0;

		//! Process a node
		virtual
		void
		emplace_back
			( ray::Node const & node
			) const = 0;

	}; // Store
#endif


#endif // Refraction_save_INCL_

