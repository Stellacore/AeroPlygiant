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


/*! \file
 *
 * \brief Unit test for class IndexVolume
 *
 */


#include "tst.hpp"

#include "env.hpp"


namespace
{
	//! Test construction of IndexVolume with no argument
	struct TestEmpty : public env::IndexVolume
	{
		inline
		TestEmpty
			()
			: IndexVolume()
		{ }

		inline
		double
		nuValue
			( engabra::g3::Vector const & rVec
			) const
		{
			return 1.;
		}

	}; // TestEmpty

	//! Test construction of IndexVolume with argument
	struct TestVolume : public env::IndexVolume
	{
		inline
		TestVolume
			()
			: IndexVolume(env::sAllSpace)
		{ }

		inline
		double
		nuValue
			( engabra::g3::Vector const & rVec
			) const
		{
			return 1.;
		}

	}; // TestVolume

	//! Check basic construction
	void
	test0
		( std::ostream & oss
		)
	{
		// [DoxyExample00]
		// [DoxyExample00]

		// Successful compile is test condition
		TestVolume const tVolume{};
		TestEmpty const tEmpty{};
	}
}


/*! \brief Unit test for IndexVolume
 */
int
main
	()
{
	std::ostringstream oss;

	test0(oss);

	return tst::finish(oss);
}

