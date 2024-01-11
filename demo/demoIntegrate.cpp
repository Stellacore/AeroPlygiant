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

#include <Engabra>

#include <iostream>
#include <vector>


namespace
{
	using namespace engabra::g3;

	// Configuration
	constexpr double sPeriod{ 10. };
	constexpr double sTauMax{ 3.*sPeriod };
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
		Vector theAcc{ null<Vector>() };
		Vector theVel{ null<Vector>() };
		Vector thePos{ null<Vector>() };

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
		return State{ nextTau, nextAcc, nextVel, nextPos };
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


/*! \brief Solve "rotating rocket" demo problem by two methods.
 *
 * Ref explantion in file header for demoIntegrate.cpp.
 */
int
main
	()
{
	// initial state
	State state0{ 0., acceleration(0.), zero<Vector>(), zero<Vector>() };

	// propagate forward with simple Euler's method
	std::vector<State> const eStates
		{ euler::nextStates(state0, sTauMax, sTauDel) };

	// display results of Euler integration
	std::cerr << "eStates.size: " << eStates.size() << std::endl;
	for (State const & eState : eStates)
	{
		if (nearInt(eState.theTau, sTauDel))
		{
			std::cout << "eState: " << eState << '\n';
		}
	}

}

