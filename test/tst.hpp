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


#ifndef Refraction_tst_INCL_
#define Refraction_tst_INCL_


/*! \file
 *
 * \brief Utilities that support testing
 *
 */


#include <Engabra>

#include <iostream>
#include <limits>
#include <sstream>
#include <string>


/*! \brief Functions and data that support project testing.
 */
namespace tst
{
	//! \brief compare generic 'got' and 'expected' values.
	template <typename Type>
	inline
	void
	checkGotExp
		( std::ostringstream & oss
		, Type const & got
		, Type const & exp
		, std::string const & tname
		, double const & tol = std::numeric_limits<double>::epsilon()
		)
	{
		if (! engabra::g3::nearlyEquals(got, exp, tol))
		{
			using engabra::g3::io::fixed;
			Type const dif{ got - exp };
			oss << "Failure of '" << tname << "' test\n";
			oss << "exp: " << exp << '\n';
			oss << "got: " << got << '\n';
			oss << "dif: " << fixed(dif, 3u, 18u) << '\n';
		}
	}

	//! CTest/CMake main program return conventions
	struct CTest
	{
		static constexpr int pass{ 0 }; //!< all tests successful
		static constexpr int fail{ 1 }; //!< one or more test failures
	};

	/*! \brief CMake/CTest compatible exit code based on (! msg.empty()).
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

	//! Convenience version that calls finish(std::string const &).
	inline
	int
	finish
		( std::ostringstream const & oss
		)
	{
		return finish(oss.str());
	}

} // [tst]

#endif // Refraction_tst_INCL_

