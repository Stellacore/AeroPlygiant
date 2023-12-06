//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


//! \file Simple ray tracing example to test results gathering/save.


#include "env.hpp"
#include "ray.hpp"
#include "save.hpp"

#include <Engabra>

#include <iostream>
#include <sstream>


namespace tst
{
	using namespace engabra::g3;

	/*! \brief Simple example of a spherical shape with constant index.
	 *
	 * NOTE: acts as a black hole - attracts ray paths into a circular orbit!
	 *
	 */
	struct Sphere : public env::IndexVolume
	{
		Vector const theCenter{ null<Vector>() };
		double const theRadius{ null<double>() };
		double const theNuCenter{ null<double>() };
		double const theNuEdge{ null<double>() };

		//! Construct a sphere in space
		inline
		explicit
		Sphere
			( Vector const & center
			, double const & radius
			, double const & nuCenter = 1.5
			, double const & nuEdge = 1.
			)
			: IndexVolume{}
			, theCenter{ center }
			, theRadius{ radius }
			, theNuCenter{ nuCenter }
			, theNuEdge{ nuEdge }
		{
		}

		//! Index of refraction value at vector location rVec
		inline
		virtual
		double
		nuValue
			( Vector const & rVec
			) const
		{
			double nu{ theNuEdge };
			// assume linear gradient
			// nu(r) = nuCenter - frac*(nuEdge-nuCenter)
			// nu(r) = nuCenter - (r/radius)*(nuEdge-nuCenter)
			// nu(r) = nuCenter - r*(1./radius)*(nuEdge-nuCenter)
			// grad(nu) = -(1./radius)*nuEdge-nuCenter;
			double const dist{ magnitude(rVec - theCenter) };
			double const frac{ dist / theRadius };
			if (frac < 1.)
			{
				nu = frac*(theNuEdge - theNuCenter) + theNuCenter;
			}
			return nu;
		}


		//! Gradient (approximation) of Index of refraction at location rVec
		inline
		virtual
		Vector
		nuGradient
			( Vector const & rVec
			) const
		{
			Vector grad{ zero<Vector>() };
			Vector const delta{ rVec - theCenter };
			double const dist{ magnitude(delta) };
			if (dist < theRadius)
			{
				Vector const gDir{ direction(delta) };
				double const gMag{ (-1./theRadius)*(theNuEdge-theNuCenter) };
				grad = gMag * gDir;
			}
			return grad;
		}



	}; // IndexVolume


	//! \brief Thick slab of constant index of refraction
	struct Slab : public env::IndexVolume
	{
		double const theBegX{ null<double>() };
		double const theEndX{ null<double>() };
		double const theNuInside{ null<double>() };
		double const theNuOutside{ null<double>() };

		//! Value constructor
		inline
		explicit
		Slab
			( double const & begX
			, double const & endX
			, double const & nuInside = 1.5
			, double const & nuOutside = 1.
			)
			: IndexVolume{}
			, theBegX{ begX }
			, theEndX{ endX }
			, theNuInside{ nuInside }
			, theNuOutside{ nuOutside }
		{
		}


		//! \brief Index of refraction value at vector location rVec.
		inline
		virtual
		double
		nuValue
			( Vector const & rVec
			) const
		{
			double nu{ theNuOutside };
			double const & valX = rVec[0];
			if ((theBegX < valX) && (valX < theEndX))
			{
				nu = theNuInside;
			}
			return nu;
		}

	}; // Slab

} // [tst]


/*! \brief Simple ray trace example with which to check save::Path
 */
int
main
	()
{
	int istat{ 1 }; // default to failure
	std::ostringstream oss;

	using namespace engabra::g3;

/*
	// construct test sphere object
	tst::Sphere const media
		( Vector{ 5., 5., 5. }
		, 2. // radius
		, 1.5 // nu center
		, 1.0 // nu edge
		);
*/

	tst::Slab const media
		( 4.   // xBeg
		, 6.   // xEnd
		, 1.5  // nu inside
		, 1.0  // nu outside
		);

	// basic construction of path
	Vector const tBeg{ 1., 1., 1. };
	Vector const rBeg{ 0., 0., 0.  };
	Vector const stopNear{ 10., 10., 10. };
	constexpr double saveStepDist{ 1./1. }; // save this often
	save::Path path(tBeg, rBeg, stopNear, saveStepDist);

	// create tracer
//	constexpr double propStepDist{ .001 }; // integration step size
constexpr double propStepDist{ 1./2. }; // integration step size
	ray::Propagator const prop{ &media, propStepDist };

	// interact with data consumer
	prop.traceNodes(tBeg, rBeg, &path);

	std::cout << "path.size: " << path.theNodes.size() << std::endl;
	for (ray::Node const & node : path.theNodes)
	{
//		std::cout << node.infoBrief() << std::endl;
	}

oss << "Failure of tangent Path test\n";
	if (oss.str().empty())
	{
		istat = 0;
	}
	else
	{
		std::cerr << oss.str() << '\n';
	}

	return istat;
}

