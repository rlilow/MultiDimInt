#include "BindMemberFunction.hpp"

#include <array>
#include <cmath>

#include <cstring>

#include <iostream>
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

template <std::size_t DimFix, std::size_t DimInt>
MultiDimInt::Integrator<DimFix, DimInt>::Integrator (const Integrand<DimFix, DimInt>& func, const Algorithm& alg, const std::string& identifier) :
	Func(func),
	Alg(alg.clone()),
	LowerBounds(),
	UpperBounds(),
	Identifier(identifier)
{
	static_assert(DimInt != 0, "MultiDimInt::Integrator Error: Number of integration variables is zero");
}

template <std::size_t DimFix, std::size_t DimInt>
template <class Class>
MultiDimInt::Integrator<DimFix, DimInt>::Integrator (const MemberIntegrandPointer<DimFix, DimInt, Class> memberFuncPointer, Class& object, const Algorithm& alg, const std::string& identifier) :
	Func(bind_member_function_to_object(memberFuncPointer, object)),
	Alg(alg.clone()),
	LowerBounds(),
	UpperBounds(),
	Identifier(identifier)
{
	static_assert(DimInt != 0, "MultiDimInt::Integrator Error: Number of integration variables is zero");
}

template <std::size_t DimFix, std::size_t DimInt>
template <class Class>
MultiDimInt::Integrator<DimFix, DimInt>::Integrator (const ConstMemberIntegrandPointer<DimFix, DimInt, Class> constMemberFuncPointer, const Class& constObject, const Algorithm& alg, const std::string& identifier) :
	Func(bind_member_function_to_object(constMemberFuncPointer, constObject)),
	Alg(alg.clone()),
	LowerBounds(),
	UpperBounds(),
	Identifier(identifier)
{
	static_assert(DimInt != 0, "MultiDimInt::Integrator Error: Number of integration variables is zero");
}

template <std::size_t DimFix, std::size_t DimInt>
MultiDimInt::Integrator<DimFix, DimInt>::Integrator (const Integrator& otherIntegrator) :
	Func(otherIntegrator.Func),
	Alg(otherIntegrator.Alg->clone()),
	LowerBounds(otherIntegrator.LowerBounds),
	UpperBounds(otherIntegrator.UpperBounds),
	Identifier(otherIntegrator.Identifier)
{
	static_assert(DimInt != 0, "MultiDimInt::Integrator Error: Number of integration variables is zero");
}

template <std::size_t DimFix, std::size_t DimInt>
bool MultiDimInt::Integrator<DimFix, DimInt>::integrate (const std::array<double, DimFix>& argsFix, double& value, double& error) const
{
	Algorithm::Result integral = Alg->run(bind_member_function_to_object(&Integrator::algorithm_internal_integrand_for_unit_hypercube, *this), DimInt, argsFix.data());	// the 'bind_member_function_to_object' is needed to pass the member function Integrator< DimFix, DimInt >::algorithm_internal_integrand_for_unit_hypercube as a MultiDimInt::Algorithm::InternalIntegrand
	
	value = integral.Value;
	error = integral.Error;
	
	if ( integral.Failed )	// if integration fails, call error handler and return 'false'
	{
		error_handler(argsFix, integral);
		
		return false;
	}
	else	// if integration succeeds, return 'true'
	{
		return true;
	}
}

template <std::size_t DimFix, std::size_t DimInt>
bool MultiDimInt::Integrator<DimFix, DimInt>::integrate (const std::array<double, DimFix>& argsFix, const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& upperBounds, double& value, double& error) const
{
	double signFlip = 1.0;
	
	for ( std::size_t i_argInt = 0; i_argInt < DimInt; ++i_argInt )	// check relation between lower and upper integration boundaries, as the integration algorithm Alg expects the lower ones to be samller than the corresponding upper ones
	{
		const double lowerBound = lowerBounds[i_argInt];
		const double upperBound = upperBounds[i_argInt];
		
		if ( lowerBound == upperBound )	// if both are equal for any integration variable, the integral is exact 0
		{
			value = 0.0;
			error = 0.0;
			
			return true;
		}
		
		if ( lowerBound < upperBound )	// if a lower boundary is smaller than the corresponding upper boundary, copy them unchanged into LowerBounds and UpperBounds
		{
			LowerBounds[i_argInt] = lowerBound;
			UpperBounds[i_argInt] = upperBound;
		}
		else	// if a lower boundary is larger than the corresponding upper boundary, swap them before copying them into LowerBounds and UpperBounds, and make up for this by a relative sign flip
		{
			LowerBounds[i_argInt] = upperBound;
			UpperBounds[i_argInt] = lowerBound;
			
			signFlip *= -1.0;
		}
	}
	
	Algorithm::Result integral = Alg->run(bind_member_function_to_object(&Integrator::algorithm_internal_integrand_for_custom_hypercube, *this), DimInt, argsFix.data());	// the 'bind_member_function_to_object' is needed to pass the member function Integrator< DimFix, DimInt >::algorithm_internal_integrand_for_custom_hypercube as a MultiDimInt::Algorithm::InternalIntegrand
	
	integral.Value *= signFlip;	// correct the value of the integral by the overall sign flip due to swapping lower and upper integration boundaries
	
	value = integral.Value;
	error = integral.Error;
	
	if ( integral.Failed )	// if integration fails, call error handler and return 'false'
	{
		error_handler(argsFix, integral);
		
		return false;
	}
	else	// if integration succeeds, return 'true'
	{
		return true;
	}
}

template <std::size_t DimFix, std::size_t DimInt>
void MultiDimInt::Integrator<DimFix, DimInt>::integrate_without_warning (const std::array<double, DimFix>& argsFix, double& value, double& error) const
{
	Algorithm::Result integral = Alg->run(bind_member_function_to_object(&Integrator::algorithm_internal_integrand_for_unit_hypercube, *this), DimInt, argsFix.data());	// the 'bind_member_function_to_object' is needed to pass the member function Integrator< DimFix, DimInt >::algorithm_internal_integrand_for_unit_hypercube as a MultiDimInt::Algorithm::InternalIntegrand
	
	value = integral.Value;
	error = integral.Error;
}

template <std::size_t DimFix, std::size_t DimInt>
void MultiDimInt::Integrator<DimFix, DimInt>::integrate_without_warning (const std::array<double, DimFix>& argsFix, const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& upperBounds, double& value, double& error) const
{
	double signFlip = 1.0;
	
	for ( std::size_t i_argInt = 0; i_argInt < DimInt; ++i_argInt )	// check relation between lower and upper integration boundaries, as the integration algorithm Alg expects the lower ones to be samller than the corresponding upper ones
	{
		const double lowerBound = lowerBounds[i_argInt];
		const double upperBound = upperBounds[i_argInt];
		
		if ( lowerBound == upperBound )	// if both are equal for any integration variable, the integral is exact 0
		{
			value = 0.0;
			error = 0.0;
		}
		
		if ( lowerBound < upperBound )	// if a lower boundary is smaller than the corresponding upper boundary, copy them unchanged into LowerBounds and UpperBounds
		{
			LowerBounds[i_argInt] = lowerBound;
			UpperBounds[i_argInt] = upperBound;
		}
		else	// if a lower boundary is larger than the corresponding upper boundary, swap them before copying them into LowerBounds and UpperBounds, and make up for this by a relative sign flip
		{
			LowerBounds[i_argInt] = upperBound;
			UpperBounds[i_argInt] = lowerBound;
			
			signFlip *= -1.0;
		}
	}
	
	Algorithm::Result integral = Alg->run(bind_member_function_to_object(&Integrator::algorithm_internal_integrand_for_custom_hypercube, *this), DimInt, argsFix.data());	// the 'bind_member_function_to_object' is needed to pass the member function Integrator< DimFix, DimInt >::algorithm_internal_integrand_for_custom_hypercube as a MultiDimInt::Algorithm::InternalIntegrand
	
	value = integral.Value *= signFlip;	// correct the value of the integral by the overall sign flip due to swapping lower and upper integration boundaries
	error = integral.Error;
}

template <std::size_t DimFix, std::size_t DimInt>
MultiDimInt::Integrator<DimFix, DimInt>& MultiDimInt::Integrator<DimFix, DimInt>::operator= (const Integrator& otherIntegrator)
{
	Func = otherIntegrator.Func;
	*Alg = *(otherIntegrator.Alg);
	UpperBounds = otherIntegrator.UpperBounds;
	LowerBounds = otherIntegrator.LowerBounds;
	Identifier = otherIntegrator.Identifier;
	
	return *this;
}

template <std::size_t DimFix, std::size_t DimInt>
MultiDimInt::Integrator<DimFix, DimInt>::~Integrator ()
{
	delete Alg;
}

template <std::size_t DimInt>
MultiDimInt::Integrator<0, DimInt>::Integrator (const IntegrandWithoutFixedArguments<DimInt>& func, const Algorithm& alg, const std::string& identifier) :
	Func(func),
	Alg(alg.clone()),
	LowerBounds(),
	UpperBounds(),
	Identifier(identifier)
{}

template <std::size_t DimInt>
template <class Class>
MultiDimInt::Integrator<0, DimInt>::Integrator (const MemberIntegrandWithoutFixedArgumentsPointer<DimInt, Class> memberFuncPointer, Class& object, const Algorithm& alg, const std::string& identifier) :
	Func(bind_member_function_to_object(memberFuncPointer, object)),
	Alg(alg.clone()),
	LowerBounds(),
	UpperBounds(),
	Identifier(identifier)
{}

template <std::size_t DimInt>
template <class Class>
MultiDimInt::Integrator<0, DimInt>::Integrator (const ConstMemberIntegrandWithoutFixedArgumentsPointer<DimInt, Class> constMemberFuncPointer, const Class& constObject, const Algorithm& alg, const std::string& identifier) :
	Func(bind_member_function_to_object(constMemberFuncPointer, constObject)),
	Alg(alg.clone()),
	LowerBounds(),
	UpperBounds(),
	Identifier(identifier)
{}

template <std::size_t DimInt>
MultiDimInt::Integrator<0, DimInt>::Integrator (const Integrator& otherIntegrator) :
	Func(otherIntegrator.Func),
	Alg(otherIntegrator.Alg->clone()),
	LowerBounds(otherIntegrator.LowerBounds),
	UpperBounds(otherIntegrator.UpperBounds),
	Identifier(otherIntegrator.Identifier)
{}

template <std::size_t DimInt>
bool MultiDimInt::Integrator<0, DimInt>::integrate (double& value, double& error) const
{
	Algorithm::Result integral = Alg->run(bind_member_function_to_object(&Integrator::algorithm_internal_integrand_for_unit_hypercube, *this), DimInt, NULL);	// the 'bind_member_function_to_object' is needed to pass the member function Integrator< 0, DimInt >::algorithm_internal_integrand_for_unit_hypercube as a MultiDimInt::Algorithm::InternalIntegrand
	
	value = integral.Value;
	error = integral.Error;
	
	if ( integral.Failed )	// if integration fails, call error handler and return 'false'
	{
		error_handler(integral);
		
		return false;
	}
	else	// if integration succeeds, return 'true'
	{
		return true;
	}
}

template <std::size_t DimInt>
bool MultiDimInt::Integrator<0, DimInt>::integrate (const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& upperBounds, double& value, double& error) const
{
	double signFlip = 1.0;
	
	for ( std::size_t i_argInt = 0; i_argInt < DimInt; ++i_argInt )	// check relation between lower and upper integration boundaries, as the integration algorithm Alg expects the lower ones to be samller than the corresponding upper ones
	{
		const double lowerBound = lowerBounds[i_argInt];
		const double upperBound = upperBounds[i_argInt];
		
		if ( lowerBound == upperBound )	// if both are equal for any integration variable, the integral is exact 0
		{
			value = 0.0;
			error = 0.0;
			
			return true;
		}
		
		if ( lowerBound < upperBound )	// if a lower boundary is smaller than the corresponding upper boundary, copy them unchanged into LowerBounds and UpperBounds
		{
			LowerBounds[i_argInt] = lowerBound;
			UpperBounds[i_argInt] = upperBound;
		}
		else	// if a lower boundary is larger than the corresponding upper boundary, swap them before copying them into LowerBounds and UpperBounds, and make up for this by a relative sign flip
		{
			LowerBounds[i_argInt] = upperBound;
			UpperBounds[i_argInt] = lowerBound;
			
			signFlip *= -1.0;
		}
	}
	
	Algorithm::Result integral = Alg->run(bind_member_function_to_object(&Integrator::algorithm_internal_integrand_for_custom_hypercube, *this), DimInt, NULL);	// the 'bind_member_function_to_object' is needed to pass the member function Integrator< 0, DimInt >::algorithm_internal_integrand_for_custom_hypercube as a MultiDimInt::Algorithm::InternalIntegrand
	
	integral.Value *= signFlip;	// correct the value of the integral by the overall sign flip due to swapping lower and upper integration boundaries
	
	value = integral.Value;
	error = integral.Error;
	
	if ( integral.Failed )	// if integration fails, call error handler and return 'false'
	{
		error_handler(integral);
		
		return false;
	}
	else	// if integration succeeds, return 'true'
	{
		return true;
	}
}

template <std::size_t DimInt>
void MultiDimInt::Integrator<0, DimInt>::integrate_without_warning (double& value, double& error) const
{
	Algorithm::Result integral = Alg->run(bind_member_function_to_object(&Integrator::algorithm_internal_integrand_for_unit_hypercube, *this), DimInt, NULL);	// the 'bind_member_function_to_object' is needed to pass the member function Integrator< 0, DimInt >::algorithm_internal_integrand_for_unit_hypercube as a MultiDimInt::Algorithm::InternalIntegrand
	
	value = integral.Value;
	error = integral.Error;
}

template <std::size_t DimInt>
void MultiDimInt::Integrator<0, DimInt>::integrate_without_warning (const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& upperBounds, double& value, double& error) const
{
	double signFlip = 1.0;
	
	for ( std::size_t i_argInt = 0; i_argInt < DimInt; ++i_argInt )	// check relation between lower and upper integration boundaries, as the integration algorithm Alg expects the lower ones to be samller than the corresponding upper ones
	{
		const double lowerBound = lowerBounds[i_argInt];
		const double upperBound = upperBounds[i_argInt];
		
		if ( lowerBound == upperBound )	// if both are equal for any integration variable, the integral is exact 0
		{
			value = 0.0;
			error = 0.0;
		}
		
		if ( lowerBound < upperBound )	// if a lower boundary is smaller than the corresponding upper boundary, copy them unchanged into LowerBounds and UpperBounds
		{
			LowerBounds[i_argInt] = lowerBound;
			UpperBounds[i_argInt] = upperBound;
		}
		else	// if a lower boundary is larger than the corresponding upper boundary, swap them before copying them into LowerBounds and UpperBounds, and make up for this by a relative sign flip
		{
			LowerBounds[i_argInt] = upperBound;
			UpperBounds[i_argInt] = lowerBound;
			
			signFlip *= -1.0;
		}
	}
	
	Algorithm::Result integral = Alg->run(bind_member_function_to_object(&Integrator::algorithm_internal_integrand_for_custom_hypercube, *this), DimInt, NULL);	// the 'bind_member_function_to_object' is needed to pass the member function Integrator< 0, DimInt >::algorithm_internal_integrand_for_custom_hypercube as a MultiDimInt::Algorithm::InternalIntegrand
	
	value = integral.Value * signFlip;	// correct the value of the integral by the overall sign flip due to swapping lower and upper integration boundaries
	error = integral.Error;
}

template <std::size_t DimInt>
MultiDimInt::Integrator<0, DimInt>& MultiDimInt::Integrator<0, DimInt>::operator= (const Integrator& otherIntegrator)
{
	Func = otherIntegrator.Func;
	*Alg = *(otherIntegrator.Alg);
	UpperBounds = otherIntegrator.UpperBounds;
	LowerBounds = otherIntegrator.LowerBounds;
	Identifier = otherIntegrator.Identifier;
	
	return *this;
}

template <std::size_t DimInt>
MultiDimInt::Integrator<0, DimInt>::~Integrator ()
{
	delete Alg;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

////////////////////////////////////////////////////////////////////////////////////////////////////
// private

template <std::size_t DimFix, std::size_t DimInt>
void MultiDimInt::Integrator<DimFix, DimInt>::error_handler (const std::array<double, DimFix>& argsFix, const Algorithm::Result& integral) const
{
	std::cout << 															std::endl
			  << " MultiDimInt::Integrator Warning: Integration failed"	 << std::endl
			  << "	-Identifier:        " << Identifier					 << std::endl 	// cout the identifier
			  << "	-Fixed argument(s): " << argsFix[0];								// cout the fixed argument(s);
	
	for ( std::size_t i_argFix = 1; i_argFix < DimFix; ++i_argFix )
	{
		std::cout << " , " << argsFix[i_argFix];
	}
	
	std::cout <<																		std::endl
			  << "	-Value:             " << integral.Value 						 << std::endl	// cout the integral value as well as the estimated absolute and relative error
			  << "	-Absolute error:    " << integral.Error 						 << std::endl
			  << "	-Relative error:    " << std::abs(integral.Error/integral.Value) << std::endl
			  << "	------------------- " << 											std::endl
			  << integral.Comment		  << 											std::endl;	// cout the comment
}

template <std::size_t DimFix, std::size_t DimInt>
double MultiDimInt::Integrator<DimFix, DimInt>::algorithm_internal_integrand_for_unit_hypercube (const double* argsFix, const double* argsInt) const
{
	std::array<double, DimFix> argsFixStd;	// copy the contents of 'argsFix' and 'argsInt' into two std::arrays, such that they can be used as arguments of 'Func'
	std::array<double, DimInt> argsIntStd;
	
	for ( std::size_t i_argFix = 0; i_argFix < DimFix; ++i_argFix )
	{
		argsFixStd[i_argFix] = argsFix[i_argFix];
	}
	
	for ( std::size_t i_argInt = 0; i_argInt < DimInt; ++i_argInt )
	{
		argsIntStd[i_argInt] = argsInt[i_argInt];
	}
	
	return Func(argsFixStd, argsIntStd);	// call the Integrand 'Func'
}

template <std::size_t DimFix, std::size_t DimInt>
double MultiDimInt::Integrator<DimFix, DimInt>::algorithm_internal_integrand_for_custom_hypercube (const double* argsFix, const double* argsInt) const
{
	std::array<double, DimFix> argsFixStd;	// copy the contents of 'argsFix' and 'argsInt' into two std::arrays, such that they can be used as arguments of 'Func'
	std::array<double, DimInt> argsIntStd;
	
	for ( std::size_t i_argFix = 0; i_argFix < DimFix; ++i_argFix )
	{
		argsFixStd[i_argFix] = argsFix[i_argFix];
	}
	
	double jacobian = 1.0;
	
	for ( std::size_t i_argInt = 0; i_argInt < DimInt; ++i_argInt )	// thereby, perform the changes of variables necessary to implement the integral boundaries 'LowerBounds' and 'UpperBounds'
	{
		const double argInt = argsInt[i_argInt];
		
		const double lowerBound = LowerBounds[i_argInt];
		const double upperBound = UpperBounds[i_argInt];
		
		if ( lowerBound > NegativeInfinity )
		{
			if ( upperBound < PositiveInfinity )	// finite integration region (lowerBound, upperBound)
			{
				argsIntStd[i_argInt] = lowerBound + (upperBound - lowerBound) * argInt;
				
				jacobian *= upperBound - lowerBound;
			}
			else	// semi-infinite integration region (lowerBound, +infinity)
			{
				argsIntStd[i_argInt] = lowerBound - 1.0 + 1.0/argInt;
				
				jacobian *= 1.0/argInt/argInt;
			}
		}
		else
		{
			if ( upperBound < PositiveInfinity )	// semi-infinite integration region (-infinity, upperBound)
			{
				argsIntStd[i_argInt] = upperBound + 1.0 - 1.0/argInt;
				
				jacobian *= 1.0/argInt/argInt;
			}
			else	// doubly-infinite integration region (-infinity, +infinity)
			{
				argsIntStd[i_argInt] = (2.0*argInt-1.0) / argInt / (1.0-argInt);
				
				jacobian *= 1.0/argInt/argInt + 1.0/(1.0-argInt)/(1.0-argInt);
			}
		}
	}
	
	return Func(argsFixStd, argsIntStd) * jacobian;	// call the Integrand 'Func'
}

template <std::size_t DimInt>
void MultiDimInt::Integrator<0, DimInt>::error_handler (const Algorithm::Result& integral) const
{
	std::cout << 																		std::endl
			  << " MultiDimInt::Integrator Warning: Integration failed"				 << std::endl
			  << "	-Identifier:        " << Identifier								 << std::endl 	// cout the identifier, the integral value as well as the estimated absolute and relative error
			  << "	-Value:             " << integral.Value 						 << std::endl
			  << "	-Absolute error:    " << integral.Error 						 << std::endl
			  << "	-Relative error:    " << std::abs(integral.Error/integral.Value) << std::endl
			  << "	------------------- " << 											std::endl
			  << integral.Comment		  << 											std::endl;	// cout the comment
}

template <std::size_t DimInt>
double MultiDimInt::Integrator<0, DimInt>::algorithm_internal_integrand_for_unit_hypercube (const double* dummyArgsFix, const double* argsInt) const
{
	std::array<double, DimInt> argsIntStd;	// copy the content of 'argsInt' into a std::array, such that it can be used as the argument of 'Func'
	
	for ( std::size_t i_argInt = 0; i_argInt < DimInt; ++i_argInt )
	{
		argsIntStd[i_argInt] = argsInt[i_argInt];
	}
	
	return Func(argsIntStd);	// call the IntegrandWithoutFixedArguments 'Func'
}

template <std::size_t DimInt>
double MultiDimInt::Integrator<0, DimInt>::algorithm_internal_integrand_for_custom_hypercube (const double* dummyArgsFix, const double* argsInt) const
{
	std::array<double, DimInt> argsIntStd;	// copy the content of 'argsInt' into a std::array, such that it can be used as the argument of 'Func'
	
	double jacobian = 1.0;
	
	for ( std::size_t i_argInt = 0; i_argInt < DimInt; ++i_argInt )	// thereby, perform the changes of variables necessary to implement the integral boundaries 'LowerBounds' and 'UpperBounds'
	{
		const double argInt = argsInt[i_argInt];
		
		const double lowerBound = LowerBounds[i_argInt];
		const double upperBound = UpperBounds[i_argInt];
		
		if ( lowerBound > NegativeInfinity )
		{
			if ( upperBound < PositiveInfinity )	// finite integration region (lowerBound, upperBound)
			{
				argsIntStd[i_argInt] = lowerBound + (upperBound - lowerBound) * argInt;
				
				jacobian *= upperBound - lowerBound;
			}
			else	// semi-infinite integration region (lowerBound, +infinity)
			{
				argsIntStd[i_argInt] = lowerBound - 1.0 + 1.0/argInt;
				
				jacobian *= 1.0/argInt/argInt;
			}
		}
		else
		{
			if ( upperBound < PositiveInfinity )	// semi-infinite integration region (-infinity, upperBound)
			{
				argsIntStd[i_argInt] = upperBound + 1.0 - 1.0/argInt;
				
				jacobian *= 1.0/argInt/argInt;
			}
			else	// doubly-infinite integration region (-infinity, +infinity)
			{
				argsIntStd[i_argInt] = (2.0*argInt-1.0) / argInt / (1.0-argInt);
				
				jacobian *= 1.0/argInt/argInt + 1.0/(1.0-argInt)/(1.0-argInt);
			}
		}
	}
	
	return Func(argsIntStd) * jacobian;	// call the IntegrandWithoutFixedArguments 'Func'
}