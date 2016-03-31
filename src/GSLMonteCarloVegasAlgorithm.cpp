#include "GSLMonteCarloVegasAlgorithm.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::GSLMonteCarloVegasAlgorithm::GSLMonteCarloVegasAlgorithm (const double absErr, const double relErr, const std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType) :
	GSLMonteCarloAlgorithm(absErr, relErr, numEval, randomNumberGeneratorType),
	Alpha(1.5),
	Iterations(5),
	Mode(1)
{}

MultiDimInt::GSLMonteCarloVegasAlgorithm::GSLMonteCarloVegasAlgorithm (const double absErr, const double relErr, const std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType,
																	   const double alpha, const std::size_t iterations, const int mode) :
	GSLMonteCarloAlgorithm(absErr, relErr, numEval, randomNumberGeneratorType),
	Alpha(alpha),
	Iterations(iterations),
	Mode(mode)
{
	if ( Alpha < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloVegasAlgorithm Error: Value of 'alpha' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Iterations <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloVegasAlgorithm Error: Value of 'iterations' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( (Mode < -1) || (Mode > 1) )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLMonteCarloVegasAlgorithm Error: Invalid value of 'mode'" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::GSLMonteCarloVegasAlgorithm* MultiDimInt::GSLMonteCarloVegasAlgorithm::clone () const
{
	return new GSLMonteCarloVegasAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

int MultiDimInt::GSLMonteCarloVegasAlgorithm::gsl_mc_integration (const std::size_t dimInt, gsl_monte_function& gslMonteIntegrand, double& value, double& error, std::string& furtherComment) const
{
	std::vector<double> lowerBound (dimInt, 0.0);	// this must not be const, as gsl_monte_vegas_integrate expects the integration boundaries as arrays of non-const double values
	std::vector<double> upperBound (dimInt, 1.0);
	
	gsl_monte_vegas_state* workspace = gsl_monte_vegas_alloc(dimInt);
	
	gsl_monte_vegas_params params;	// set Vegas specific parameters
	
	gsl_monte_vegas_params_get(workspace, &params);
	
	params.alpha = Alpha;
	params.iterations = Iterations;
	params.mode = Mode;
	
	gsl_monte_vegas_params_set(workspace, &params);
	
	int fail = gsl_monte_vegas_integrate(&gslMonteIntegrand, lowerBound.data(), upperBound.data(), dimInt, NumEval, RandomNumberGenerator, workspace, &value, &error);
	
	gsl_monte_vegas_free(workspace);
	
	if ( (error / value > RelErr) && (error > AbsErr) && (fail == 0) )	// set fail to 14 (GSL_ETOL) if the required tolerance was not reached, but only if there has been no other error
	{
		fail = 14;
		
		const double chisq = gsl_monte_vegas_chisq(workspace);	// chi-squared per degree of freedom
	
		std::stringstream chisqComment;	// turn chisq into a string with 2-digit mantissa
		chisqComment.precision(2);
		chisqComment << std::fixed << chisq;
			
		furtherComment = std::string("\n")
					   + std::string("	-Chi-squared per degree of freedom: ") + chisqComment.str();
	}
	
	return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private