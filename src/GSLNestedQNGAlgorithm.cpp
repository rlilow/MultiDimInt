#include "GSLNestedQNGAlgorithm.h"

#include <gsl/gsl_integration.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::GSLNestedQNGAlgorithm::GSLNestedQNGAlgorithm (const double absErr, const double relErr) :
	GSLNestedAlgorithm(absErr, relErr)
{}

MultiDimInt::GSLNestedQNGAlgorithm* MultiDimInt::GSLNestedQNGAlgorithm::clone () const
{
	return new GSLNestedQNGAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

int MultiDimInt::GSLNestedQNGAlgorithm::gsl_integration (const gsl_function& recursiveIntegrand, double& value, double& error) const
{
	std::size_t neval;	// number of evaluations needed by the GSL QNG routine (will not be used)
	
	const int fail = gsl_integration_qng(&recursiveIntegrand, 0.0, 1.0, AbsErr, RelErr, &value, &error, &neval);
	
	return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private