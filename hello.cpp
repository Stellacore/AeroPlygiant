//
// Copyright (c) 2023 Stellacore Corporation. All rights reserved.
//


#include <Engabra>

#include <iostream>

int
main
	()
{
	std::cout << "Hello Refraction\n";
	std::cout << "Using Engabra Project: " << engabra::projectVersion() << '\n';
	std::cout << "Using Engabra Source: " << engabra::sourceIdentity() << '\n';
}
