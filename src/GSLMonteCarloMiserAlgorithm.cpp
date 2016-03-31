#include "GSLMonteCarloMiserAlgorithm.h"

#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_miser.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::GSLMonteCarloMiserAlgorithm::GSLMonteCarloMiserAlgorithm (const double absErr, const double relErr, const std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType) :
	GSLMonteCarloAlgorithm(absErr, relErr, numEval, randomNumberGeneratorType),
	Estimate_frac(0.1),
	Min_calls_per_dim(16),
	Min_calls_per_dim_per_bisection(512),
	Alpha(2.0),
	Dither(0.0)
{}

MultiDimInt::GSLMonteCarloMiserAlgorithm::GSLMonteCarloMiserAlgorithm (const double absErr, const double relErr, const std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType,
																	   const double estimate_frac, const std::size_t min_calls_per_dim, const std::size_t min_calls_per_dim_per_bisection, const double alpha, const double dither) :
	GSLMonteCarloAlgorithm(absErr, relErr, numEval, randomNumberGeneratorType),
	Estimate_frac(estimate_frac),
	Min_calls_per_dim(min_calls_per_dim),
	Min_calls_per_dim_per_bisection(min_calls_per_dim_per_bisection),
	Alpha(alpha),
	Dither(dither)
{
	if ( Estimate_frac <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloMiserAlgorithm Error: Value of 'estimate_frac' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	
	if ( Min_calls_per_dim <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloMiserAlgorithm Error: Value of 'min_calls_per_dim' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Min_calls_per_dim_per_bisection <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloMiserAlgorithm Error: Value of 'min_calls_per_dim_per_bisection' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Alpha <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloMiserAlgorithm Error: Value of 'alpha' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Dither < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloMiserAlgorithm Error: Value of 'dither' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::GSLMonteCarloMiserAlgorithm* MultiDimInt::GSLMonteCarloMiserAlgorithm::clone () const
{
	return new GSLMonteCarloMiserAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

int MultiDimInt::GSLMonteCarloMiserAlgorithm::gsl_mc_integration (const std::size_t dimInt, gsl_monte_function& gslMonteIntegrand, double& value, double& error, std::string& furtherComment) const
{
	std::vector<double> lowerBound (dimInt, 0.0);	// this must not be const, as gsl_monte_miser_integrate expects the integration boundaries as arrays of non-const double values
	std::vector<double> upperBound (dimInt, 1.0);
	
	gsl_monte_miser_state* workspace = gsl_monte_miser_alloc(dimInt);
	
	gsl_monte_miser_params params;	// set Miser specific parameters
	
	gsl_monte_miser_params_get(workspace, &params);
	
	params.estimate_frac = Estimate_frac;
	params.min_calls = Min_calls_per_dim * dimInt;
	params.min_calls_per_bisection = Min_calls_per_dim_per_bisection * dimInt;
	params.alpha = Alpha;
	params.dither = Dither;
	
	gsl_monte_miser_params_set(workspace, &params);
	
	int fail = gsl_monte_miser_integrate(&gslMonteIntegrand, lowerBound.data(), upperBound.data(), dimInt, NumEval, RandomNumberGenerator, workspace, &value, &error);
	
	gsl_monte_miser_free(workspace);
	
	if ( (error / value > RelErr) && (error > AbsErr) && (fail == 0) )	// set fail to 14 (GSL_ETOL) if the required tolerance was not reached, but only if there has been no other error
	{
		fail = 14;
	}
	
	return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private