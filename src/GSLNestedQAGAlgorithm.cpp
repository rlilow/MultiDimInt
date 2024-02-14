#include "GSLNestedQAGAlgorithm.hpp"

#include <iostream>

#include <gsl/gsl_integration.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::GSLNestedQAGAlgorithm::GSLNestedQAGAlgorithm (const double absErr, const double relErr, const std::size_t maxInterval, const int key) :
	GSLNestedAlgorithm(absErr, relErr),
	MaxInterval(maxInterval),
	Key(key)
{
	if ( MaxInterval <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLNestedQAGAlgorithm Error: Maximal number of intervals is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( (Key < 1) || (Key > 6) )
	{
		std::cout << std::endl
				  << " MultiDimInt::GSLNestedQAGAlgorithm Error: Key out of range" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::GSLNestedQAGAlgorithm* MultiDimInt::GSLNestedQAGAlgorithm::clone () const
{
	return new GSLNestedQAGAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

int MultiDimInt::GSLNestedQAGAlgorithm::gsl_integration (const gsl_function& recursiveIntegrand, double& value, double& error) const
{
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(MaxInterval);
	
	const int fail = gsl_integration_qag(&recursiveIntegrand, 0.0, 1.0, AbsErr, RelErr, MaxInterval, Key, workspace, &value, &error);
	
	gsl_integration_workspace_free(workspace);
	
	return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private