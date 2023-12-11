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

#ifndef aply_ray_PathView_INCL_
#define aply_ray_PathView_INCL_

/*! \file
 *
 * \brief Ray propagation simulation functions.
 *
 */


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

	//! Provide view of interesting path information
	struct PathView
	{
		//! Must be set by consumer
		std::vector<Node> const * const thePtNodes;
		std::vector<double> const * const thePtArcDists;

		/*! \brief Attach an instance to (EXTERNALLY managed!!) path data.
		 *
		 * \note Does *NOT* take position of the path data, but only
		 * accesses the externally owned and managed instance.
		 */
		explicit
		PathView
			( Path const * const & ptPath
			)
			: thePtNodes{ &(ptPath->theNodes) }
			, thePtArcDists{ &(ptPath->theArcDists) }
		{ }

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

		//! Distance along path (propagation resolution approximation).
		inline
		double
		pathDistance
			() const
		{
			double const distSum
				{ std::accumulate
					(thePtArcDists->cbegin(), thePtArcDists->cend(), 0.)
				};
			return distSum;
		}

		//! Distance subtended by begDeviation() at pathDistance().
		inline
		double
		begDeflection
			() const
		{
			return (magnitude(begDeviation()) * pathDistance());
		}

		//! Distance subtended by endDeviation() at pathDistance().
		inline
		double
		endDeflection
			() const
		{
			return (magnitude(endDeviation()) * pathDistance());
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
			oss << '\n';
			oss << "  pathDistance: " << pathDistance();
			oss << '\n';
			oss << " begDeflection: " << io::fixed(begDeflection(), 3u, 3u);
			oss << '\n';
			oss << " endDeflection: " << io::fixed(endDeflection(), 3u, 3u);
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
} // [aply]


#endif // aply_ray_PathView_INCL_

