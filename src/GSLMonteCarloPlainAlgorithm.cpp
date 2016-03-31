#include "GSLMonteCarloPlainAlgorithm.h"

#include <string>
#include <vector>

#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_rng.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::GSLMonteCarloPlainAlgorithm::GSLMonteCarloPlainAlgorithm (const double absErr, const double relErr, const std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType) :
	GSLMonteCarloAlgorithm(absErr, relErr, numEval, randomNumberGeneratorType)
{}

MultiDimInt::GSLMonteCarloPlainAlgorithm* MultiDimInt::GSLMonteCarloPlainAlgorithm::clone () const
{
	return new GSLMonteCarloPlainAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

int MultiDimInt::GSLMonteCarloPlainAlgorithm::gsl_mc_integration (const std::size_t dimInt, gsl_monte_function& gslMonteIntegrand, double& value, double& error, std::string& furtherComment) const
{
	std::vector<double> lowerBound (dimInt, 0.0);	// this must not be const, as gsl_monte_plain_integrate expects the integration boundaries as arrays of non-const double values
	std::vector<double> upperBound (dimInt, 1.0);
	
	gsl_monte_plain_state* workspace = gsl_monte_plain_alloc(dimInt);
	
	int fail = gsl_monte_plain_integrate(&gslMonteIntegrand, lowerBound.data(), upperBound.data(), dimInt, NumEval, RandomNumberGenerator, workspace, &value, &error);
	
	gsl_monte_plain_free(workspace);
	
	if ( (error / value > RelErr) && (error > AbsErr) && (fail == 0) )	// set fail to 14 (GSL_ETOL) if the required tolerance was not reached, but only if there has been no other error
	{
		fail = 14;
	}
	
	return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private