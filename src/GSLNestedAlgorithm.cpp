#include "GSLNestedAlgorithm.hpp"

#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::Algorithm::Result MultiDimInt::GSLNestedAlgorithm::run (const InternalIntegrand& func, const std::size_t dimInt, const double* argsFix) const
{
	Algorithm::Result integral = {false, 0.0, 0.0, ""};
	
	double* argsInt = new double[dimInt];
	
	GSLNestedData gslNestedData(0, func, argsFix, argsInt, dimInt, this); // '0' is the initital value of 'NestingCounter'
	
	gsl_set_error_handler_off();	// turn off the GSL error handler during the integration, as Integrator::error_handler takes care of errors
	
	gsl_function recursiveIntegrand = {&internal_recursion, &gslNestedData};	// integrating over the outermost integration variable
	
	int fail;	// if the integration succeeds, this is set to '0', otherwise it is a non-vanishing GSL error code
	
	fail = gsl_integration(recursiveIntegrand, integral.Value, integral.Error);
	
	gsl_set_error_handler(NULL);	// turn on the default GSL error handler again
	
	delete[] argsInt;
	
	if ( fail != 0 )	// if integration failed, write the GSL error message and 'furtherComment' into 'integral.Comment'
	{
		integral.Failed = true;
		
		integral.Comment = std::string("	-GSL error: ") + std::string(gsl_strerror(fail));
	}
	
	return integral;
}

bool MultiDimInt::GSLNestedAlgorithm::is_parallelized () const
{
  return false;	// all nested GSL integration algorithms use only a single core
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

MultiDimInt::GSLNestedAlgorithm::GSLNestedAlgorithm (const double absErr, const double relErr) :
	Algorithm(absErr, relErr)
{}

double MultiDimInt::GSLNestedAlgorithm::internal_recursion (const double currentArg, void* gslNestedData)
{
	GSLNestedData data = *(GSLNestedData*) gslNestedData;	// it is important to create a new instance instead of just passing the pointer, as only this allows to keep track of the correct current recursion depth
	
	data.ArgsInt[data.NestingCounter] = currentArg;	// inserting the current one-dimensional integration variable 'currentArg' into the respective element of ArgsInt
	
	data.NestingCounter++;
	
	double result;
	
	if ( data.NestingCounter < data.DimInt )
	{
		gsl_function recursiveIntegrand = {&internal_recursion, &data};	// integrating over the next integration variable (one recursion level deeper)
		
		double error;	// error estimated by the GSL routine for inner integrations gets discarded
		
		data.ThisGSLNestedAlgorithm->gsl_integration(recursiveIntegrand, result, error);
	}
	else
	{
		result = data.Func(data.ArgsFix, data.ArgsInt);	// at the innermost recursion level the integrand gets evaluated at the current integration variables stored in 'ArgsInt'
	}
	
	return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private