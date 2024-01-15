#include "CubaAlgorithm.h"

#include <bitset>
#include <climits>
#include <iostream>
#include <sstream>
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::Algorithm::Result MultiDimInt::CubaAlgorithm::run (const InternalIntegrand& func, const std::size_t dimInt, const double* argsFix) const
{
	if ( dimInt > INT_MAX )	// Cuba algorithms only accept an 'int' as the number of integration variables, not a potentially larger 'size_t'; in practice 'dimInt' will of course never exceed the largest possible 'int' value, but explicitly checking this won't hurt
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaAlgorithm Error: Number of integration variables to large" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	Algorithm::Result integral = {false, 0.0, 0.0, ""};
	
	CubaData cubaData(func, argsFix);
	
	bool integrationSucceeded;		// if the integration succeeds, this is set to 'true', otherwise to 'false'
	double prob;					// chi^2 probability that the estimated error is not a reliable estimate of the true integration error
	std::string furtherComment("");	// if the integration fails, an additional comment may be written to this string
	
	integrationSucceeded = cuba_integration(static_cast<int>(dimInt), cubaData, integral.Value, integral.Error, prob, furtherComment);
	
	if ( not integrationSucceeded )	// if integration failed, write Cuba error message, the value of 'prob', and 'furtherComment' into 'integral.Comment'
	{
		integral.Failed = true;
		
		std::stringstream probComment;	// turn 'prob' into a string with 2-digit mantissa
		probComment.precision(2);
		probComment << std::fixed << prob;
		
		integral.Comment = std::string("	-Cuba error: Failed to reach the specified tolerance") + std::string("\n")
						 + std::string("	-Probability that the error estimate is not reliable: ") + probComment.str()
						 + furtherComment;
	}
	
	return integral;
}

bool MultiDimInt::CubaAlgorithm::is_parallelized () const
{
  return true;	// all Cuba integration algorithms parallelize the sampling of the integrand
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

MultiDimInt::CubaAlgorithm::CubaAlgorithm (const double absErr, const double relErr, const int maxEval) :
	Algorithm(absErr, relErr),
	MaxEval(maxEval),
	Flags(0),
	Mineval(0),
	Statefile(),
	Spin(NULL)
{
	if ( MaxEval <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaAlgorithm Error: Maximal number of integrand evaluations is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::CubaAlgorithm::CubaAlgorithm (const double absErr, const double relErr, const int maxEval,
										   const int flags, const int mineval, const std::string& statefile, void* spin) :
	Algorithm(absErr, relErr),
	MaxEval(maxEval),
	Flags(flags),
	Mineval(mineval),
	Statefile(statefile),
	Spin(spin)
{
	if ( MaxEval <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaAlgorithm Error: Maximal number of integrand evaluations is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( (Flags < 0) )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaAlgorithm Error: Invalid value of 'flags'" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Mineval < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaAlgorithm Error: Value of 'mineval' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

int MultiDimInt::CubaAlgorithm::cuba_integrand (const int* dimInt, const double* argsInt, const int* ncomp, double* result, void* cubaData)
{
	CubaData* data = (CubaData*) cubaData;
	
	result[0] = data->Func(data->ArgsFix, argsInt);	// evaluate integrand
	
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private