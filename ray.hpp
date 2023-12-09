//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
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
		static
		Start
		from
			( Vector const & anyTan
			, Vector const & loc
			)
		{
			return { direction(anyTan), loc };
		}

	}; // Start


	//! Characterization of ray path tangent interacting at step boundary.
	enum DirChange
	{
		  Null      //!< Unset
		, Unaltered //!< Tangent dir unchanged (no gradient)
		, Converged //!< Tangent dir refracted toward gradient (into dense)
		, Diverged  //!< Tangent dir refracted away from gradient (into sparse)
		, Reflected //!< Tangent dir reflected from boundary (total internal)

	}; // DirChange

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
		return revChange;
	}


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
		infoBrief
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
		infoString
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
		reversed
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
		, env::IndexVolume const * const & ptMedia //!< Refraction media
		)
	{
		// default to unaltered
		std::pair<Vector, DirChange> tanDirChange{ tDirPrev, Unaltered };
		Vector & tDirNext = tanDirChange.first;
		DirChange & tChange = tanDirChange.second;
		//
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
			if (nuPrev < nuNext) // propagating into more dense media
			{
				Spinor const spin{  rootXi, currB };
				tDirNext = (spin * gCurrInv).theVec;
				tChange = Converged;
			}
			else
			if (nuNext < nuPrev) // propagating into less dense media
			{
				Spinor const spin{ -rootXi, currB };
				tDirNext = (spin * gCurrInv).theVec;
				tChange = Diverged;
			}
			// (nuNext == nuPrev) // same as default (gCurr == 0)
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

		struct Step
		{
			double theNextNu;
			Vector theNextTan;
			DirChange theChange;

		}; // Step

		//! Estimate next tangent based on local object refraction
		inline
		Step
		nextStep
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
std::ostringstream oss;
oss << "gCurr: " << gCurr;
			if (! (gTol < gMag)) // unaltered
			{
oss << " PassThrough";
				// location at which to evaluate nuNext (iteratively refined)
				Vector const qNext{ rCurr + .5*theStepDist*tNext };

				// update estimated forward next refraction index value
				// (which should be same as previous)
				nuNext = thePtMedia->nuValue(qNext);

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
//std::cout << "tNext: " << io::fixed(tNext) << std::endl;
//std::cout << "gCurr: " << io::fixed(gCurr) << std::endl;
					if (! isReflection)
					{
oss << " Refraction ";
						// update refraction index to midpoint of predicted next
						// interval (along evolving next tangent direction).
						qNext = rCurr + .5*theStepDist*tNext;
					}
					else
					// if (isReflection) // perfect reflection
					{
oss << " REFLECTION ";
						Vector const gDir{ direction(gCurr) };
						qNext = rCurr + .5*theStepDist*gDir;
						// No need to iterate further for perfect reflection
						doLoop = false; // exit after updating next ray step
					}

					// update estimated forward next refraction index value
					nuNext = thePtMedia->nuValue(qNext);

					// update tangent direction
					std::pair<Vector, DirChange> const tDirChange
						{ nextTangentDir
							(tPrev, nuPrev, gCurr, nuNext, thePtMedia)
						};
					isReflection = (Reflected == tDirChange.second);

					// evaluate convergence of tangent direction
					Vector const & tResult = tDirChange.first;
					difSq = magSq(tResult - tNext);

					// update tangent direction
					tNext = tResult;
					change = tDirChange.second;

				} // while loop on refraction index estimation

/*
if (0u < numLoop)
{
std::cout
	<< "-- numLoop: " << numLoop
	<< "  at: " << rCurr
	<< '\n';
}
*/

			} // significant gCurr magnitude

oss << " nuPrev: " << io::fixed(nuPrev, 3u, 3u);
oss << " nuNext: " << io::fixed(nuNext, 3u, 3u);
oss << " tNext: " << tNext;
//std::cout << oss.str() << std::endl;

			return Step{ nuNext, tNext, change };
		}

		//! Predicted next location stepsize units along tangent from rVec
		inline
		Vector
		nextLocation
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
		isValid
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
		traceNodes
			( Start const & start
			, Consumer * const & ptConsumer
			) const
		{
			if (isValid())
			{
				// start with initial conditions
				Vector tPrev{ start.theTanDir };
				Vector rCurr{ start.thePntLoc };

				// incident media IoR
				Vector const & tBeg = start.theTanDir;
				Vector const & rBeg = start.thePntLoc;
				double nuPrev
					{ thePtMedia->nuValue(rBeg - .5*theStepDist*tBeg) };

				// propagate until path approximate reaches requested length
				while (ptConsumer->size() < ptConsumer->capacity())
				{
					// determine propagation change at this step
					Step const stepNext{ nextStep(tPrev, nuPrev, rCurr) };
					Vector const & tNext = stepNext.theNextTan;
					double const & nuNext = stepNext.theNextNu;
					DirChange const & change = stepNext.theChange;

					// propagate ray to next node location
					Vector const rNext{ nextLocation(rCurr, tNext) };

					// give consumer opportunity to record node data
					ptConsumer->emplace_back
						(Node{ tPrev, nuPrev, rCurr, nuNext, tNext, change });

					// update state for next node
					tPrev = tNext;
					rCurr = rNext;
					nuPrev = nuNext;
				}
//std::cout << "\n####\n#### done tracing\n####\n";
			}
		}

/*TODO - unused*/
		//! \brief Return nodes in a std::vector structure.
		inline
		std::vector<Node>
		nodePath
			( Start const & start
			, double const & nominalLength
			) const
		{
			std::vector<Node> nodes;

			if (isValid())
			{
				// estimate return structure size and allocate space
				constexpr double padFactor{ 9./8. }; // about 12% extra
				double const dubSize{ padFactor * nominalLength / theStepDist };
				std::size_t const nomSize{ static_cast<std::size_t>(dubSize) };

				std::size_t stride{ 1u };
				std::size_t saveSize{ nomSize };
				if (1000u < nomSize)
				{
					stride = nomSize / 1000u;
				}
				std::size_t const useSize{ nomSize / stride + 10u };
				nodes.reserve(useSize);

				traceNodes(start, &nodes);
std::cout << "nodes.size; " << nodes.size() << std::endl;
			}

			return nodes;
		}
/*TODO*/

	}; // Propagator


	/*! \brief Store path information suitable for later visualization.
	 *
	 * Basic methods are compatible with those of std::vector.
	 *
	 */
	struct Path
	{
		Start const theStart{};
		Vector const theStopLoc{ null<Vector>() };
		double const theSaveDelta{ null<double>() };
		std::vector<ray::Node> theNodes{};
		double thePrevNearDist{ null<double>() };
		double theCurrNearDist{ null<double>() };

		//! Construct storage based on nominal distance between points
		inline
		explicit
		Path
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
			, theCurrNearDist{ magnitude(theStopLoc - theStart.thePntLoc) }
		{
			// estimate distance (as if straight line)
			double const nomDist{ magnitude(stopNearTo - theStart.thePntLoc) };
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
				addNode(node);
			}
		}

	}; // Path

} // [ray]

#endif // Refraction_ray_INCL_

