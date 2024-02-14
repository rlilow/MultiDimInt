#include "GSLMonteCarloAlgorithm.hpp"

#include <iostream>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::Algorithm::Result MultiDimInt::GSLMonteCarloAlgorithm::run (const InternalIntegrand& func, const std::size_t dimInt, const double* argsFix) const
{
	Algorithm::Result integral{false, 0.0, 0.0, ""};
	
	GSLMonteCarloData gslMonteCarloData(func, argsFix);
	
	int fail;						// if the integration succeeds, this is set to 0, otherwise it contains a GSL error code
	std::string furtherComment("");	// if the integration fails, an additional comment may be written to this string
	
	gsl_monte_function gslMonteIntegrand = {&gsl_mc_integrand, dimInt, &gslMonteCarloData};
	
	fail = gsl_mc_integration(dimInt, gslMonteIntegrand, integral.Value, integral.Error, furtherComment);
	
	if ( fail != 0 )	// if integration failed, write the GSL error message and 'furtherComment' into 'integral.Comment'
	{
		integral.Failed = true;
		
		integral.Comment = std::string("	-GSL error: ") + std::string(gsl_strerror(fail))
						 + furtherComment;
	}
	
	return integral;
}

bool MultiDimInt::GSLMonteCarloAlgorithm::is_parallelized () const
{
  return false;	// all GSL Monte Carlo integration algorithms use only a single core
}

MultiDimInt::GSLMonteCarloAlgorithm& MultiDimInt::GSLMonteCarloAlgorithm::operator= (const GSLMonteCarloAlgorithm& otherGSLMonteCarloAlgorithm)
{
	Algorithm::operator=(otherGSLMonteCarloAlgorithm);	// calling the assignment operator of the base class
	
	NumEval = otherGSLMonteCarloAlgorithm.NumEval;
	RandomNumberGeneratorType = otherGSLMonteCarloAlgorithm.RandomNumberGeneratorType;
	
	gsl_rng_free (RandomNumberGenerator);								// free the old GSL random number generator...
	RandomNumberGenerator = gsl_rng_alloc(RandomNumberGeneratorType);	// ...and allocate it with the new type
	
	gsl_rng_set(RandomNumberGenerator, 0);
	
	return *this;
}

MultiDimInt::GSLMonteCarloAlgorithm::~GSLMonteCarloAlgorithm ()
{
	gsl_rng_free (RandomNumberGenerator);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

MultiDimInt::GSLMonteCarloAlgorithm::GSLMonteCarloAlgorithm (const double absErr, const double relErr, const std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType) :
	Algorithm(absErr, relErr),
	NumEval(numEval),
	RandomNumberGeneratorType(randomNumberGeneratorType)
{
	if ( NumEval <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloAlgorithm Error: Number of integrand evaluations is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	RandomNumberGenerator = gsl_rng_alloc(RandomNumberGeneratorType);
	
	gsl_rng_set(RandomNumberGenerator, 0);
}

MultiDimInt::GSLMonteCarloAlgorithm::GSLMonteCarloAlgorithm (const GSLMonteCarloAlgorithm& otherGSLMonteCarloAlgorithm) :
	Algorithm(otherGSLMonteCarloAlgorithm),
	NumEval(otherGSLMonteCarloAlgorithm.NumEval),
	RandomNumberGeneratorType(otherGSLMonteCarloAlgorithm.RandomNumberGeneratorType)
{
	RandomNumberGenerator = gsl_rng_alloc(RandomNumberGeneratorType);
	
	gsl_rng_set(RandomNumberGenerator, 0);
}

double MultiDimInt::GSLMonteCarloAlgorithm::gsl_mc_integrand (double* argsInt, const size_t dimInt, void* gslMonteCarloData)
{
	GSLMonteCarloData* data = (GSLMonteCarloData*) gslMonteCarloData;
	
	return data->Func(data->ArgsFix, argsInt);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private