#include "Algorithm.hpp"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

MultiDimInt::Algorithm::Algorithm (const double absErr, const double relErr) :
	AbsErr(absErr),
	RelErr(relErr)
{
	if ( AbsErr < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::Algorithm Error: Absolute error limit is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( RelErr < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::Algorithm Error: Relative error limit is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private