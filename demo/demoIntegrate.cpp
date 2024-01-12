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

/*! \file \brief Demonstrate numerical integration for simple ODE.
 *
 * Demonstrate use of math::DiffEqSystem and math::DiffEqSolve to compute
 * the solution of a simple second order differential equation system.
 *
 * The problem is associated with a object undergoing constant acceleration
 * for which the magnitude is constant, but for which the direction is
 * changing uniformly over time. (Ref demoIntegrate.lyx)
 *
 * The code in this file solves the problem with two approaches. The first
 * is a simple fininte difference estimation (essentially Euler's method).
 * The second solution uses the 4th order Runge Kutta algorithm from 
 * math::DiffEqSolve (with a vector valued system of equations).
 *
*/

#include "mathDiffEqSolve.hpp"
#include "mathDiffEqSystem.hpp"

#include <Engabra>

#include <iostream>
#include <vector>


namespace
{
	using namespace engabra::g3;

	// Configuration
	constexpr double sPeriod{ 8. };
	constexpr double sTauMax{ 4.*sPeriod };
	constexpr double sTauDel{ sPeriod/1024. };

	//! Rotating acceleration vector
	inline
	Vector
	acceleration
		( double const & tau
		)
	{
		Vector acc{ zero<Vector>() };
		// if (! (sPeriod < tau))
		{
			static double const omegaMag{ 2.*pi/sPeriod };
			static BiVector const omegaDir{ e12 };
			acc = (e1 * exp(tau * omegaMag * omegaDir)).theVec;
		}
		return acc;
	}

	//! Motion state (time, position, velocity, acceleration)
	struct State
	{
		double theTau{ null<double>() };
		Vector thePos{ null<Vector>() };
		Vector theVel{ null<Vector>() };
		Vector theAcc{ null<Vector>() };

	}; // State

	//! Standard stream output of motion State
	inline
	std::ostream &
	operator<<
		( std::ostream & strm
		, State const & state
		)
	{
		strm
			<< "tau: " << io::fixed(state.theTau, 2u, 9u)
			<< "  acc: " << io::fixed(state.theAcc, 2u, 9u)
			<< "  vel: " << io::fixed(state.theVel, 2u, 9u)
			<< "  pos: " << io::fixed(state.thePos, 2u, 9u)
			;
		return strm;
	}

	//! Difference (component by component) between two states
	inline
	State
	operator-
		( State const & stateA
		, State const & stateB
		)
	{
		return State
			{ stateA.theTau - stateB.theTau
			, stateA.thePos - stateB.thePos
			, stateA.theVel - stateB.theVel
			, stateA.theAcc - stateB.theAcc
			};
	}

	//! True if value is closest to (positive side of) an integer
	bool
	nearInt
		( double const & value
		, double const & valDelta
		)
	{
		bool isNear{ false };
		double const currFrac{ value - std::floor(value) };
		double const prevFrac{ value - std::floor(value - valDelta) };
		isNear = (currFrac < prevFrac);
		return isNear;
	}

} // [anon]

/* \brief Finite step solution (Euler's forward method 1st order).
 */
namespace euler
{

	//! Apply finite forward update
	inline
	State
	nextState
		( State const & curr
		, double const & dtau
		)
	{
		double const & tau = curr.theTau;
		Vector const currAcc{ acceleration(tau) };
		Vector const nextVel{ curr.theVel + currAcc * dtau };
		Vector const nextPos{ curr.thePos + nextVel * dtau };
		double const nextTau{ tau + dtau };
		Vector const nextAcc{ acceleration(nextTau) };
		return State{ nextTau, nextPos, nextVel, nextAcc };
	}


	//! Propagate states forward
	inline
	std::vector<State>
	nextStates
		( State const & state0 //!< Initial state
		, double const & tauMax //!< Propagate until this time
		, double const & tauDel //!< Step size (in time)
		)
	{
		std::vector<State> states; // trajectory to be computed

		std::size_t const numSamps
			{ static_cast<std::size_t>(tauMax / tauDel) + 1u };
		states.reserve(numSamps);


		// start with initial state
		states.emplace_back(state0);

		// integrate forward step by step
		for (std::size_t nSamp{1u} ; nSamp < numSamps ; ++nSamp)
		{
			State const next{ nextState(states.back(), tauDel) };
			states.emplace_back(next);
		}

		return states;
	}

} // [euler]

/*! \brief Runge-Kutta (RK) solution approach.
 */
namespace rk
{
	//! \brief System of vector equations for solving acceleration() ODE.
	struct AccelSystem : public aply::math::DiffEqSystem
	{
		double const theInitTau{ null<double>() };
		Vector const theInitPos{ null<Vector>() };
		Vector const theInitVel{ null<Vector>() };

		explicit
		AccelSystem
			( double const initTau
			, Vector const initPos
			, Vector const initVel
			)
			: DiffEqSystem()
			, theInitTau{ initTau }
			, theInitPos{ initPos }
			, theInitVel{ initVel }
		{ }

		virtual ~AccelSystem() = default;

		/*! \brief TODO
		*
		* Functions are:
		* \arg y0c1 position, e1 component
		* \arg y0c2 position, e2 component
		* \arg y0c3 position, e3 component
		* \arg y1c1 velocity, e1 component
		* \arg y1c2 velocity, e2 component
		* \arg y1c3 velocity, e3 component
		* \arg y2c1 acceleration, e1 component (from accel model)
		* \arg y2c2 acceleration, e2 component (from accel model)
		* \arg y2c3 acceleration, e3 component (from accel model)
		*
		* Derivative are:
		* y1 = y0'
		* y2 = y1'
 		*/
		virtual
		std::vector<double>
		operator()
			( std::pair<double, std::vector<double> > const & tyValues
			) const
		{
			double const & tau = tyValues.first;
			std::vector<double> const & yFuncs = tyValues.second;
			double const * ptFuncComp = yFuncs.data();
			double const & y0c1 = *ptFuncComp++;
			double const & y0c2 = *ptFuncComp++;
			double const & y0c3 = *ptFuncComp++;
			double const & y1c1 = *ptFuncComp++;
			double const & y1c2 = *ptFuncComp++;
			double const & y1c3 = *ptFuncComp++;

			Vector const acc{ acceleration(tau) };

			double const y0c1Prime = y1c1;
			double const y0c2Prime = y1c2;
			double const y0c3Prime = y1c3;
			double const y1c1Prime = acc[0];
			double const y1c2Prime = acc[1];
			double const y1c3Prime = acc[2];

			return
				{ // y0Prime
				  y0c1Prime
				, y0c2Prime
				, y0c3Prime
				  // y1Prime
				, y1c1Prime
				, y1c2Prime
				, y1c3Prime
				};
		}

		/*! \brief Initial conditios for spinning rocket problem
 		 *
 		 * Init Conditions
		 * \arg Position (y0c[012] = 0.)
		 * \arg Velocity (y1c[012] = 0.)
 		 */
		std::pair<double, std::vector<double> >
		initValues
			() const
		{
			return
				{ theInitTau
				, std::vector<double>
					{  // Pos(t0)
					  theInitPos[0]
					, theInitPos[1]
					, theInitPos[2]
					   // Vel(t0)
					, theInitVel[0]
					, theInitVel[1]
					, theInitVel[2]
					}
				};
		}

	}; // AccelSystem

	//! Use RK4 solver to approximation solution for one display step.
	inline
	State
	nextState
		( State const & currState //!< Initial state
		, double const & nextTau //!< Propagate until this time
		, double const & tauDel //!< Step size (in time)
		)
	{
		// integrate until the next reporting step

		// setup system at start of this step
		AccelSystem const accelSystem
			(currState.theTau, currState.thePos, currState.theVel);

		// solve system until next step
		aply::math::DiffEqSolve solver(tauDel);

		std::pair<double, std::vector<double> >
			const soln{ solver.solutionFor(nextTau, accelSystem) };

		std::vector<double> const & sVals = soln.second;
		Vector const nextPos{ sVals[0], sVals[1], sVals[2] };
		Vector const nextVel{ sVals[3], sVals[4], sVals[5] };

		return State{ nextTau, nextPos, nextVel, acceleration(nextTau) };
	}

	//! Use RK4 solver to approximation solution until tauMax
	inline
	std::vector<State>
	nextStates
		( State const & state0 //!< Initial state
		, double const & tauMax //!< Propagate until this time
		, double const & tauDel //!< Step size (in time)
		)
	{
		std::vector<State> states;

		// set initial state
		State currState{ state0 };
		states.emplace_back(currState);

		// iterate over future states
		// (solving each state with a many step RK4 integration)
		while (! (tauMax < currState.theTau))
		{
			// integrate until the next reporting step
			double const nextTau{ currState.theTau + tauDel };

			// record this way point
			currState = nextState(currState, nextTau, tauDel);
			states.emplace_back(currState);
		}

		return states;
	}

} // [rk]


/*! \brief Solve "rotating rocket" demo problem by two methods.
 *
 * Ref explantion in file header for demoIntegrate.cpp.
 */
int
main
	()
{
	// initial state
	State state0{ 0., zero<Vector>(), zero<Vector>(), acceleration(0.) };

	// propagate forward with simple Euler's method
	std::vector<State> const eStates
		{ euler::nextStates(state0, sTauMax, sTauDel) };

	// approximate with 4th order RK solution.
	std::vector<State> const rkStates
		{ rk::nextStates(state0, sTauMax, sTauDel) };

	// display results of Euler integration
	std::cerr << "\neStates.size: " << eStates.size() << std::endl;
	for (State const & eState : eStates)
	{
		if (nearInt(eState.theTau, sTauDel))
		{
			std::cout << "eState: " << eState << '\n';
		}
	}

	// display results of Euler integration
	std::cerr << "\nrkStates.size: " << rkStates.size() << std::endl;
	for (State const & rkState : rkStates)
	{
		if (nearInt(rkState.theTau, sTauDel))
		{
			std::cout << "rkState: " << rkState << '\n';
		}
	}

	std::size_t const eSize{ eStates.size() };
	std::size_t const rkSize{ rkStates.size() };
	std::size_t const minSize{ std::min(eSize, rkSize) };
	for (std::size_t nn{0u} ; nn < minSize ; nn += 128u)
	{
		State difState{ eStates[nn] - rkStates[nn] };
		difState.theTau = eStates[nn].theTau;
		std::cout << "difState: " << difState << '\n';
	}

}

