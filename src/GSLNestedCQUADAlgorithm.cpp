#include "GSLNestedCQUADAlgorithm.h"

#include <iostream>

#include <gsl/gsl_integration.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::GSLNestedCQUADAlgorithm::GSLNestedCQUADAlgorithm (const double absErr, const double relErr, const std::size_t maxInterval) :
	GSLNestedAlgorithm(absErr, relErr),
	MaxInterval(maxInterval)
{
	if ( MaxInterval <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLNestedCQUADAlgorithm Error: Maximal number of stored intervals is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::GSLNestedCQUADAlgorithm* MultiDimInt::GSLNestedCQUADAlgorithm::clone () const
{
	return new GSLNestedCQUADAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

int MultiDimInt::GSLNestedCQUADAlgorithm::gsl_integration (const gsl_function& recursiveIntegrand, double& value, double& error) const
{
	gsl_integration_cquad_workspace* workspace = gsl_integration_cquad_workspace_alloc(MaxInterval);
	
	const int fail = gsl_integration_cquad(&recursiveIntegrand, 0.0, 1.0, AbsErr, RelErr, workspace, &value, &error, NULL);
	
	gsl_integration_cquad_workspace_free(workspace);
	
	return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private