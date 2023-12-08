//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


//! \file Utilities that support testing


#include <iostream>
#include <sstream>
#include <string>


namespace tst
{
	//! CTest/CMake main program return conventions
	struct CTest
	{
		static constexpr int pass{ 0 }; //!< all tests successful
		static constexpr int fail{ 1 }; //!< one or more test failures
	};

	/*! CMake/CTest compatible exit code based on condition (! msg.empty())
	 *
	 * If message string is empty, return success.
	 *
	 * Example Usage Convention:
	 * \snippet test_test.cpp DoxyExample00
	 */
	inline
	int
	finish
		( std::string const & msg
		)
	{
		int istat{ CTest::fail };
		if (msg.empty())
		{
			istat = CTest::pass;
		}
		else
		{
			std::cerr << msg << '\n';
		}
		return istat;
	}

	//! Convenience version that calls finish(std::string const &)
	inline
	int
	finish
		( std::ostringstream const & oss
		)
	{
		return finish(oss.str());
	}

} // [tst]

