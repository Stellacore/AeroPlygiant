

#include <Engabra>

#include <iostream>
#include <vector>


namespace
{
	using namespace engabra::g3;

	//
	// Configuration
	//
	constexpr double sPeriod{ 10. };
	constexpr double sTauMax{ 3.*sPeriod };
	constexpr double sTauDel{ sPeriod/1024. };

	struct State
	{
		double theTau{ null<double>() };
		Vector theAcc{ null<Vector>() };
		Vector theVel{ null<Vector>() };
		Vector thePos{ null<Vector>() };

	}; // State

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

	//! Rotating acceleration vector
	inline
	Vector
	acceleration
		( double const & tau
		)
	{
		Vector acc{ zero<Vector>() };
	//	if (! (sPeriod < tau))
		{
			static double const omegaMag{ 2.*pi/sPeriod };
			static BiVector const omegaDir{ e12 };
			acc = (e1 * exp(tau * omegaMag * omegaDir)).theVec;
		}
		return acc;
	}

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

}


int
main
	()
{
	std::vector<State> states;

	std::size_t const numSamps
		{ static_cast<std::size_t>(sTauMax / sTauDel) + 1u };
	states.reserve(numSamps);

	State prev{ 0., acceleration(0.), zero<Vector>(), zero<Vector>() };

	std::cout << "prev: " << prev << '\n';
	states.emplace_back(prev);
	for (std::size_t nSamp{1u} ; nSamp < numSamps ; ++nSamp)
	{
		State const next{ nextState(states.back(), sTauDel) };

		std::cout << "next: " << next << '\n';
		states.emplace_back(next);
	}

	std::cerr << "size: " << states.size() << std::endl;

}

