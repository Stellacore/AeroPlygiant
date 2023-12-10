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

#ifndef Refraction_ray_INCL_
#define Refraction_ray_INCL_

/*! \file
 *
 * \brief Ray propagation simulation functions.
 *
 */


#include "env.hpp"

#include <Engabra>

#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


/*! \brief Functions and classes for simulation of ray propagation
 */
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


	//! Update tangent direction across single (idealized) interface boundary.
	inline
	std::pair<Vector, DirChange>
	nextTangentDir
		( Vector const & tDirPrev //!< Must be unit length
		, double const & nuPrev //!< Incoming IoR
		, Vector const & gCurr //! Must be non-zero (to be invertable
		, double const & nuNext //!< Exiting IoR
		)
	{
		// default case is an unaltered ray
		std::pair<Vector, DirChange> tanDirChange{ tDirPrev, Unaltered };
		Vector & tDirNext = tanDirChange.first;
		DirChange & tChange = tanDirChange.second;
		//
		// check for stop condition
		if (! engabra::g3::isValid(nuPrev))
		{
			tChange = Stopped;
		}
		else
		{
			// compute refraction bivector
			// note that magnitude is order of |gCurr|
			BiVector const currB{ (nuPrev/nuNext) * (tDirPrev*gCurr).theBiv };
			//
			// note that sq(bivector) = -magSq(bivector)
			double const gCurrSq{ magSq(gCurr) };
			double const radicand{ gCurrSq - magSq(currB) };
			//
			// use current conditions to select computation option
			//
			Vector const gCurrInv{ (1./gCurrSq) * gCurr };
			if (radicand < 0.) // total internal reflection
			{
				// reflect tangent from interface plane (dual to gCurr)
				tDirNext = -(gCurr * tDirPrev * gCurrInv).theVec;
				tChange = Reflected;
			}
			else
			{
				double const rootXi{ std::sqrt(radicand) };
				double const tDotG{ (tDirPrev * gCurr).theSca[0] };
				if (tDotG < 0.) // propagating into less dense media
				{
					Spinor const spin{ -rootXi, currB };
					tDirNext = (spin * gCurrInv).theVec;
					tChange = Diverged;
				}
				else
				if (0. < tDotG) // propagating into more dense media
				{
					Spinor const spin{  rootXi, currB };
					tDirNext = (spin * gCurrInv).theVec;
					tChange = Converged;
				}
				// (nuNext == nuPrev) // same as default (gCurr == 0)
			}
		}
		//
		return tanDirChange;
	}

	//! Ray propagation functions
	struct Propagator
	{
		env::IndexVolume const * const thePtMedia{ nullptr };
		double const theStepDist{ null<double>() };

	private:

		struct Step // Propagator::
		{
			double theNextNu;
			Vector theNextTan;
			DirChange theChange;

		}; // Step

		//! Estimate next tangent based on local object refraction
		inline
		Step
		nextStep // Propagator::
			( Vector const & tPrev //!< Must be unit length
			, double const & nuPrev
			, Vector const & rCurr
			) const
		{
			double nuNext{ null<double>() };
			Vector tNext{ tPrev }; // iteratively evolved from here
			DirChange change{ Null };
			//
			// Check if there's anything to compute (vs unaltered propagation)
			Vector const gCurr{ thePtMedia->nuGradient(rCurr, theStepDist) };
			double const gMag{ magnitude(gCurr) };
			static double const gTol // enough to unitize and invert gCurr
				{ std::numeric_limits<double>::min() };
			if (! (gTol < gMag)) // unaltered
			{
				// location at which to evaluate nuNext (iteratively refined)
				Vector const qNext{ rCurr + .5*theStepDist*tNext };

				// update estimated forward next refraction index value
				// (which should be same as previous)
				nuNext = thePtMedia->qualifiedNuValue(qNext);

				// update tangent direction
				tNext = tPrev; // default initialized

				// update boundary action
				change = Unaltered; // default initialized
			}
			else
			// if (gTol < gMag) // ray path (tangent direction) changes
			{
				bool isReflection{ false }; // used to exit early 
				//
				// iterate on determination of exit media IoR
				// (ref Refraction.lyx doc)
				//
				// location at which to evaluate nuNext (iteratively refined)
				Vector qNext{ rCurr + .5*theStepDist*tNext };
				//
				// tolerance until epsilon < difSq (sqrt(eps)<|dif|)
				static double const tolDifSq
					{ std::numeric_limits<double>::epsilon() };
				double difSq{ 2.*tolDifSq }; // large value to force first loop
				constexpr std::size_t maxLoop{ 10u }; // avoid infinite loop
				std::size_t numLoop{ 0u };
				bool doLoop{ true };
				while ((tolDifSq < difSq) && (numLoop++ < maxLoop) && doLoop)
				{
					// update IoR evaluation location
					if (! isReflection)
					{
						// update refraction index to midpoint of predicted next
						// interval (along evolving next tangent direction).
						qNext = rCurr + .5*theStepDist*tNext;
					}
					else
					// if (isReflection) // perfect reflection
					{
						Vector const gDir{ direction(gCurr) };
						qNext = rCurr + .5*theStepDist*gDir;
						// No need to iterate further for perfect reflection
						doLoop = false; // exit after updating next ray step
					}

					// update estimated forward next refraction index value
					nuNext = thePtMedia->qualifiedNuValue(qNext);

					// update tangent direction
					std::pair<Vector, DirChange> const tDirChange
						{ nextTangentDir(tPrev, nuPrev, gCurr, nuNext) };

					// check for stop condition
					if (Stopped == tDirChange.second)
					{
						break;
					}

					// note reflection condition for next loop iteration
					isReflection = (Reflected == tDirChange.second);

					// evaluate convergence of tangent direction
					Vector const & tResult = tDirChange.first;
					difSq = magSq(tResult - tNext);

					// update tangent direction
					tNext = tResult;
					change = tDirChange.second;

				} // while loop on refraction index estimation

			} // significant gCurr magnitude

			// check for invalid media index volume (e.g. exit region)
			if (! engabra::g3::isValid(nuNext))
			{
				change = Stopped;
			}

			return Step{ nuNext, tNext, change };
		}

		//! Predicted next location stepsize units along tangent from rVec
		inline
		Vector
		nextLocation // Propagator::
			( Vector const & rVec
			, Vector const & tVec
			) const
		{
			return { rVec + theStepDist * tVec };
		}

	public:

		//! True if this instance is valid
		inline
		bool
		isValid // Propagator::
			() const
		{
			return engabra::g3::isValid(theStepDist);
		}

		/*! Perform forward integration step by step.
		 *
		 * Essentially is Euler's method for integration of the ray path
		 * (with all attendent pitfalls).
		 *
		 */
		template <typename Consumer>
		inline
		void
		tracePath // Propagator::
			( Consumer * const & ptConsumer
			) const
		{
			if (isValid() && ptConsumer)
			{
				Vector const & tBeg = ptConsumer->theStart.theTanDir;
				Vector const & rBeg = ptConsumer->theStart.thePntLoc;

				// start with initial conditions
				Vector tPrev{ tBeg };
				Vector rCurr{ rBeg };

				// incident media IoR
				Vector const rPrev{ (rBeg - .5*theStepDist*tBeg) };
				double nuPrev{ thePtMedia->qualifiedNuValue(rPrev) };

				// propagate until path approximate reaches requested length
				// or encounteres a NaN value for index of refraction
				while (ptConsumer->size() < ptConsumer->capacity())
				{
					// determine propagation change at this step
					Step const stepNext{ nextStep(tPrev, nuPrev, rCurr) };
					Vector const & tNext = stepNext.theNextTan;
					double const & nuNext = stepNext.theNextNu;
					DirChange const & change = stepNext.theChange;

					// check for ray termination condition
					if (Stopped == change)
					{
						break;
					}

					// propagate ray to next node location
					Vector const rNext{ nextLocation(rCurr, tNext) };

					// give consumer opportunity to record node data
					Node const nextNode
						{ tPrev, nuPrev, rCurr, nuNext, tNext, change };
					ptConsumer->emplace_back(nextNode);

					// update state for next node
					tPrev = tNext;
					rCurr = rNext;
					nuPrev = nuNext;
				}
			}
		}

	}; // Propagator


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

	//! Provide view of interesting path information
	struct PathView
	{
		//! Must be set by consumer
		std::vector<Node> const * const thePtNodes;

		//! First node in path (or null if not present)
		inline
		Node
		begNode // PathView::
			() const
		{
			if (! thePtNodes->empty())
			{
				return thePtNodes->front();
			}
			return Node{};
		}

		//! Last node in path (or null if not present)
		inline
		Node
		endNode // PathView::
			() const
		{
			if (! thePtNodes->empty())
			{
				return thePtNodes->back();
			}
			return Node{};
		}

		//! Direction (of tangent) at first node
		inline
		Vector
		begDirection // PathView::
			() const
		{
			return begNode().thePrevTan;
		}

		//! Direction (of tangent) at last node
		inline
		Vector
		endDirection // PathView::
			() const
		{
			return endNode().theNextTan;
		}

		//! Direction of direct path (from first location to end location)
		inline
		Vector
		netDirection // PathView::
			() const
		{
			Vector dir{ null<Vector>() };
			if (1u < thePtNodes->size())
			{
				Vector const netDiff
					{ endNode().theCurrLoc - begNode().theCurrLoc };
				dir = direction(netDiff);
			}
			
				return dir;
		}

		//! Directed angle between fromVec and intoVec
		inline
		BiVector
		angleFromInto // PathView::
			( Vector const & fromVec
			, Vector const & intoVec
			) const
		{
			BiVector angle{ null<BiVector>() };
			if (engabra::g3::isValid(fromVec) && engabra::g3::isValid(intoVec))
			{
				Spinor const expSpin(fromVec * intoVec);
				Spinor const logSpin{ logG2(expSpin) };
				angle = logSpin.theBiv;
			}
			return angle;
		}

		//! Angle from netDirection() toward begin tangent
		inline
		BiVector
		begDeviation // PathView::
			() const
		{
			return angleFromInto(netDirection(), begDirection());
		}

		//! Angle from netDirection() toward end tangent
		inline
		BiVector
		endDeviation // PathView::
			() const
		{
			return angleFromInto(netDirection(), endDirection());
		}

		//! Angle from begDir toward endDir
		inline
		BiVector
		totalDeviation // PathView::
			() const
		{
			return angleFromInto(begDirection(), endDirection());
		}

		//! Summary of overall path info
		inline
		std::string
		infoCurvature // PathView::
			() const
		{
			std::ostringstream oss;
			oss << "begDirection: " << begDirection();
			oss << '\n';
			oss << "endDirection: " << endDirection();
			oss << '\n';
			oss << "  begDeviation: " << begDeviation();
			oss << '\n';
			oss << "  endDeviation: " << endDeviation();
			oss << '\n';
			oss << "totalDeviation: " << totalDeviation();
			return oss.str();
		}

		//! Summary of overall path info
		inline
		std::string
		infoShape // PathView::
			() const
		{
			std::ostringstream oss;
			oss << "begNode: " << begNode().infoString();
			oss << '\n';
			oss << "endNode: " << endNode().infoString();
			oss << '\n';
			oss << infoCurvature() << '\n';
			return oss.str();
		}

	}; // PathView

} // [ray]


namespace
{

	//! Overload output for ray::Start
	inline
	std::ostream &
	operator<<
		( std::ostream & ostrm
		, ray::Start const & start
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
		, ray::Node const & node
		)
	{
		ostrm << node.infoString();
		return ostrm;
	}

} // [anon]


#endif // Refraction_ray_INCL_

