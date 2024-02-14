#include "GSLNestedQAGSAlgorithm.hpp"

#include <iostream>

#include <gsl/gsl_integration.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::GSLNestedQAGSAlgorithm::GSLNestedQAGSAlgorithm (const double absErr, const double relErr, const std::size_t maxInterval) :
	GSLNestedAlgorithm(absErr, relErr),
	MaxInterval(maxInterval)
{
	if ( MaxInterval <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLNestedQAGSAlgorithm Error: Maximal number of intervals is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::GSLNestedQAGSAlgorithm* MultiDimInt::GSLNestedQAGSAlgorithm::clone () const
{
	return new GSLNestedQAGSAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

int MultiDimInt::GSLNestedQAGSAlgorithm::gsl_integration (const gsl_function& recursiveIntegrand, double& value, double& error) const
{
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(MaxInterval);
	
	const int fail = gsl_integration_qags(&recursiveIntegrand, 0.0, 1.0, AbsErr, RelErr, MaxInterval, workspace, &value, &error);
	
	gsl_integration_workspace_free(workspace);
	
	return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private