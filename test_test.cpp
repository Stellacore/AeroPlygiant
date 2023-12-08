//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


//! \file Example unit testing template


#include "tst.hpp"


// [DoxyExample00]
/*! \brief Example of test program structure/conventions.
 */
int
main
	()
{
	std::ostringstream oss;

	bool const testIsGood{ oss.good() }; // dummy for example
	if (! testIsGood)
	{
		oss << "Failure information goes here\n";
		oss << "(in whatever form with whatever diagnostics)\n";
	}
	// else // if (! testIsGood)
	// {
	// 	// Passing tests do *NOT* put anything in the error message string.
	// }

	return tst::finish(oss);
}
// [DoxyExample00]

