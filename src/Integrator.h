#ifndef MULTIDIMINT_INTEGRATOR_H
#define MULTIDIMINT_INTEGRATOR_H

#include "Algorithm.h"

#include <array>
#include <functional>
#include <limits>
#include <string>

namespace MultiDimInt
{
	/**
	 * The arguments of integrands are expected to be specified as a \c double \c std::array of length \a Dim.
	 */
	template <std::size_t Dim>
	using Arguments = std::array<double, Dim>;
	
	/**
	 * Functions that only shall be integrated over some of their arguments are expected to be of this form: They take a
	 * reference to some MultiDimInt::Arguments \a argsFix containing the fixed arguments of the function and a reference
	 * to some MultiDimInt::Arguments \a argsInt containing the integration variables, and return a \c double.
	 * 
	 * Note that the use of \c std::function allows to pass even function objects as well as lambda functions to an
	 * Integrator as long as they are of this form.
	 */
	template <std::size_t DimFix, std::size_t DimInt>
	using Integrand = std::function<double(const Arguments<DimFix>& argsFix, const Arguments<DimInt>& argsInt)>;
	
	/**
	 * Pointers to member functions that only shall be integrated over some of their arguments are expected to be of
	 * this form: They point to a method of some class \a Class which takes a reference to some MultiDimInt::Arguments \a argsFix
	 * containing the fixed arguments of the function and a reference to some MultiDimInt::Arguments \a argsInt containing
	 * the integration variables, and return a \c double.
	 * 
	 * Note that one additionally needs an instance of the host class \a Class to pass a member function to an Integrator.
	 */
	template <std::size_t DimFix, std::size_t DimInt, class Class>
	using MemberIntegrandPointer = double(Class::*)(const Arguments<DimFix>& argsFix, const Arguments<DimInt>& argsInt);
	
	/**
	 * Pointers to \c const member functions that only shall be integrated over some of their arguments are expected
	 * to be of this form: They point to a \c const method of some class \a Class which takes a reference to some
	 * MultiDimInt::Arguments \a argsFix containing the fixed arguments of the function and a reference to some MultiDimInt::Arguments
	 * \a argsInt containing the integration variables, and returns a \c double.
	 */
	template <std::size_t DimFix, std::size_t DimInt, class Class>
	using ConstMemberIntegrandPointer = double(Class::*)(const Arguments<DimFix>& argsFix, const Arguments<DimInt>& argsInt) const;
	
	/**
	 * Functions that shall be integrated over all of their arguments are expected to be of this form: They take a reference
	 * to a MultiDimInt::Arguments \a argsInt containing all the variables, and return a \c double.
	 * 
	 * Note that the use of \c std::function allows to pass even function objects as well as lambda functions to an Integrator
	 * as long as they are of this form.
	 */
	template <std::size_t DimInt>
	using IntegrandWithoutFixedArguments = std::function<double(const Arguments<DimInt>& argsInt)>;
	
	/**
	 * Pointers to member functions that shall be integrated over all of their arguments are expected to be of this
	 * form: They point to a method of some class \a Class which takes some MultiDimInt::Arguments<\a DimInt> containing
	 * all the variables, and returns a \c double.
	 */
	template <std::size_t DimInt, class Class>
	using MemberIntegrandWithoutFixedArgumentsPointer = double(Class::*)(const Arguments<DimInt>& argsInt);
	
	/**
	 * Pointers to \c const member functions that shall be integrated over all of their arguments are expected to be
	 * of this form: They point to a \c const method of some class \a Class which takes some MultiDimInt::Arguments \a argsInt
	 * containing all the variables, and returns a \c double.
	 */
	template <std::size_t DimInt, class Class>
	using ConstMemberIntegrandWithoutFixedArgumentsPointer = double(Class::*)(const Arguments<DimInt>& argsInt) const;
	
	/**
	 * Use this to specify a positive infinite integration boundary.
	 */
	constexpr double PositiveInfinity = std::numeric_limits<double>::max();
	
	/**
	 * Use this to specify a negative infinite integration boundary.
	 */
	constexpr double NegativeInfinity = std::numeric_limits<double>::lowest();
	
	/**
	 * \brief Class performing a multi-dimensional integration over some arguments of a given function.
	 *
	 * It uses an Algorithm of choice and can deal with different finite as well as infinite integration boundaries.
	 * 
	 * It expects the number of fixed arguments \a DimFix as well as the number of the integration variables \a DimInt
	 * of the integrand as template parameters.
	 */
	template <std::size_t DimFix, std::size_t DimInt>
	class Integrator
	{
	public:
		/**
		 * Constructor instantiating an integrator that integrates over the \a \DimInt integration variables of the
		 * MultiDimInt::Integrand \a func while keeping its other \a DimFix arguments fixed, using the integration Algorithm
		 * \a alg. Furthermore, one can provide an optional \c string \a identifier that will be used by Integrator::error_handler
		 * and can be useful to distinguish the warning messages of several Integrator objects used in parallel.
		 */
		Integrator (const Integrand<DimFix, DimInt>& func, const Algorithm& alg, const std::string& identifier = "");
		
		/**
		 * Constructor instantiating an integrator that integrates over the \a \DimInt integration variables of the
		 * member function characterized by the MultiDimInt::MemberIntegrandPointer \a memberFuncPointer as well as
		 * the instance \a object of its host class \a Class, while keeping its other \a DimFix arguments fixed, using
		 * the integration Algorithm \a alg. Furthermore, one can provide an optional \c string \a identifier that will
		 * be used by Integrator::error_handler and can be useful to distinguish the warning messages of several Integrator
		 * objects used in parallel.
		 */
		template <class Class>
		Integrator (const MemberIntegrandPointer<DimFix, DimInt, Class> memberFuncPointer, Class& object, const Algorithm& alg, const std::string& identifier = "");
		
		/**
		 * Constructor instantiating an integrator that integrates over the \a \DimInt integration variables of the
		 * \c const member function characterized by the MultiDimInt::ConstMemberIntegrandPointer \a constMemberFuncPointer
		 * as well as the \c const instance \a constObject of its host class \a Class, while keeping its other \a DimFix
		 * arguments fixed, using the integration Algorithm \a alg. Furthermore, one can provide an optional \c string
		 * \a identifier that will be used by Integrator::error_handler and can be useful to distinguish the warning messages
		 * of several Integrator objects used in parallel.
		 */
		template <class Class>
		Integrator (const ConstMemberIntegrandPointer<DimFix, DimInt, Class> constMemberFuncPointer, const Class& constObject, const Algorithm& alg, const std::string& identifier = "");
		
		/**
		 * Copy-constructor taking care of properly copying the integration Algorithm pointed to by Integrator::Alg
		 * from the Integrator \a otherIntegrator.
		 */
		Integrator (const Integrator& otherIntegrator);
		
		/**
		 * Writes the result of performig the integral for fixed arguments \a argsFix over the unit hypercube, i.e. integrating
		 * each variable from 0 to 1, into \a value and its estimated absolute error into \a error. If the integration succeeds,
		 * it returns \c true. Otherwise it returns \c false and also calls Integrator::error_handler.
		 */
		bool integrate (const Arguments<DimFix>& argsFix, double& value, double& error) const;
		
		/**
		 * Writes the result of performig the integral for fixed arguments \a argsFix over the hypercube specified by \a lowerBounds
		 * and \a uppperBounds into \a value and its estimated absolute error into \a error. Infinite integration boundaries
		 * can be specified using MultiDimInt::NegativeInfinity and MultiDimInt::PositiveInfinity. If the integration succeeds,
		 * it returns \c true. Otherwise it returns \c false and also calls Integrator::error_handler.
		 * 
		 * Note that this is slightly slower than using Integrator::integrate(double& value, double& error) const, which
		 * integrates over the unit hypercube. For maximal performance we thus recommend using that method instead and taking
		 * care of the integration boundaries by applying appropriate changes of variables in the integrand Integrator::Func.
		 */
		bool integrate (const Arguments<DimFix>& argsFix, const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const;
		
		/**
		 * Writes the result of performig the integral for fixed arguments \a argsFix over the unit hypercube, i.e. integrating
		 * each variable from 0 to 1, into \a value and its estimated absolute error into \a error. In contrast to
		 * Integrator::integrate(double& value, double& error), this method does not check if the integration succeeded
		 * or failed.
		 */
		void integrate_without_warning (const Arguments<DimFix>& argsFix, double& value, double& error) const;
		
		/**
		 * Writes the result of performig the integral for fixed arguments \a argsFix over the hypercube specified by \a lowerBounds
		 * and \a uppperBounds into \a value and its estimated absolute error into \a error. Infinite integration boundaries
		 * can be specified using MultiDimInt::NegativeInfinity and MultiDimInt::PositiveInfinity. In contrast to
		 * Integrator::integrate(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error),
		 * this method does not check if the integration succeeded or failed.
		 * 
		 * Note that this is slightly slower than using Integrator::integrate_without_warning(double& value, double& error),
		 * which integrates over the unit hypercube. For maximal performance we thus recommend using that method instead
		 * and taking care of the integration boundaries by applying appropriate changes of variables in the integrand Integrator::Func.
		 */
		void integrate_without_warning (const Arguments<DimFix>& argsFix, const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const;
		
		/**
		 * Assignment operator taking care of properly copying the integration Algorithm pointed to by Integrator::Alg from
		 * the Integrator \a otherIntegrator.
		 */
		Integrator& operator= (const Integrator& otherIntegrator);
		
		/**
		 * Destructor deleting the integration Algorithm pointed to by Integrator::Alg.
		 */
		~Integrator ();
		
	private:
		/**
		 * Function that shall be integrated.
		 */
		Integrand<DimFix, DimInt> Func;
		
		/**
		 * Pointer to the integration Algorithm.
		 */
		Algorithm* Alg;
		
		/**
		 * Lower integration boundaries.
		 * 
		 * These have to be \c mutable, as they will be modified by the methods
		 * Integrator::integrate(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const and
		 * Integrator::integrate_without_warning(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const.
		 */
		mutable Arguments<DimInt> LowerBounds;
		
		/**
		 * Upper integration boundaries.
		 * 
		 * These have to be \c mutable, as they will be modified by the methods
		 * Integrator::integrate(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const and
		 * Integrator::integrate_without_warning(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const.
		 */
		mutable Arguments<DimInt> UpperBounds;
		
		/**
		 * The \c string that will be written to the standard output by Integrator::error_handler to identify the specific
		 * Integrator object.
		 */
		std::string Identifier;
		
		/**
		 * Is called when the integration run performed by Integrator::integrate(double& value, double& error) const or
		 * Integrator::integrate(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const
		 * fails. It writes Integrator::Identifier and all relevant information that can be extracted from the Algorithm::Result
		 * \a integral, i.e. the integral value, its estimated absolute and relative errors as well as possible further
		 * comments passed by the Algorithm pointed to by Integrator::Alg, to the standard output.
		 */
		void error_handler (const Arguments<DimFix>& argsFix, const Algorithm::Result& integral) const;
		
		/**
		 * Wrapper for Integrator::Func that provides the form of the integrand expected by an integration Algorithm, i.e.
		 * an Algorithm::InternalIntegrand, if the integral shall be performed over the unit hypercube.
		 */
		double algorithm_internal_integrand_for_unit_hypercube (const double* dummyArgsFix, const double* argsInt) const;
		
		/**
		 * Wrapper for Integrator::Func that provides the form of the integrand expected by an integration Algorithm, i.e.
		 * an Algorithm::InternalIntegrand, if the integral shall be performed over the hypercube specified by \a LowerBounds
		 * and \a UpperBounds.
		 */
		double algorithm_internal_integrand_for_custom_hypercube (const double* dummyArgsFix, const double* argsInt) const;
	};
	
	/**
	 * \brief Class performing a multi-dimensional integration over all arguments of a given function.
	 * 
	 * This is a specialization of Integrator, which generally allows to keep some of the arguments fixed.
	 *
	 * It uses an Algorithm of choice and can deal with different finite as well as infinite integration boundaries.
	 * 
	 * It expects the the number of variables \a DimInt of the integrand as a template parameter.
	 */
	template <std::size_t DimInt>
	class Integrator<0, DimInt>
	{
	public:
		/**
		 * Constructor instantiating an integrator that integrates over all \a \DimInt variables of the
		 * MultiDimInt::IntegrandWithoutFixedArguments \a func, using the integration Algorithm \a alg. Furthermore,
		 * one can provide an optional \c string \a identifier that will be used by Integrator<0, DimInt>::error_handler
		 * and can be useful to distinguish the warning messages of several Integrator<0, DimInt> objects used in parallel.
		 */
		Integrator (const IntegrandWithoutFixedArguments<DimInt>& func, const Algorithm& alg, const std::string& identifier = "");
		
		/**
		 * Constructor instantiating an integrator that integrates over all \a \DimInt variables of the member function
		 * characterized by the MultiDimInt::MemberIntegrandWithoutFixedArgumentsPointer \a memberFuncPointer as well as
		 * the instance \a object of its host class \a Class, using the integration Algorithm \a alg. Furthermore, one
		 * can provide an optional \c string \a identifier that will be used by Integrator<0, DimInt>::error_handler and
		 * can be useful to distinguish the warning messages of several Integrator<0, DimInt> objects used in parallel.
		 */
		template <class Class>
		Integrator (const MemberIntegrandWithoutFixedArgumentsPointer<DimInt, Class> memberFuncPointer, Class& object, const Algorithm& alg, const std::string& identifier = "");
		
		/**
		 * Constructor instantiating an integrator that integrates over the \a \DimInt variables of the \c const member
		 * function characterized by the MultiDimInt::ConstMemberIntegrandWithoutFixedArgumentsPointer \a constMemberFuncPointer
		 * as well as the \c const instance \a constObject of its host class \a Class, using the integration Algorithm
		 * \a alg. Furthermore, one can provide an optional \c string \a identifier that will be used by Integrator<0, DimInt>::error_handler
		 * and can be useful to distinguish the warning messages of several Integrator<0, DimInt> objects used in parallel.
		 */
		template <class Class>
		Integrator (const ConstMemberIntegrandWithoutFixedArgumentsPointer<DimInt, Class> constMemberFuncPointer, const Class& constObject, const Algorithm& alg, const std::string& identifier = "");
		
		/**
		 * Copy-constructor taking care of properly copying the integration Algorithm pointed to by Integrator<0, DimInt>::Alg
		 * from the Integrator<0, DimInt> \a otherIntegrator.
		 */
		Integrator (const Integrator& otherIntegrator);
		
		/**
		 * Writes the result of performig the integral over the unit hypercube, i.e. integrating each variable from 0 to
		 * 1, into \a value and its estimated absolute error into \a error. If the integration succeeds, it returns \c true.
		 * Otherwise it returns \c false and also calls Integrator<0, DimInt>::error_handler.
		 */
		bool integrate (double& value, double& error) const;
		
		/**
		 * Writes the result of performig the integral over the hypercube specified by \a lowerBounds and \a uppperBounds
		 * into \a value and its estimated absolute error into \a error. Infinite integration boundaries can be specified
		 * using MultiDimInt::NegativeInfinity and MultiDimInt::PositiveInfinity. If the integration succeeds, it returns
		 * \c true. Otherwise it returns \c false and also calls Integrator<0, DimInt>error_handler.
		 * 
		 * Note that this is slightly slower than using Integrator<0, DimInt>::integrate(double& value, double& error) const,
		 * which integrates over the unit hypercube. For maximal performance we thus recommend using that method instead
		 * and taking care of the integration boundaries by applying appropriate changes of variables in the integrand
		 * Integrator<0, DimInt>::Func.
		 */
		bool integrate (const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const;
		
		/**
		 * Writes the result of performig the integral over the unit hypercube, i.e. integrating each variable from 0 to 1,
		 * into \a value and its estimated absolute error into \a error. In contrast to Integrator<0, DimInt>::integrate(double& value, double& error) const,
		 * this method does not check if the integration succeeded or failed.
		 */
		void integrate_without_warning (double& value, double& error) const;
		
		/**
		 * Writes the result of performig the integral over the hypercube specified by \a lowerBounds and \a uppperBounds
		 * into \a value and its estimated absolute error into \a error. Infinite integration boundaries can be specified
		 * using MultiDimInt::NegativeInfinity anf MultiDimInt::PositiveInfinity. In contrast to
		 * Integrator<0, DimInt>::integrate(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const,
		 * this method does not check if the integration succeeded or failed.
		 * 
		 * Note that this is slightly slower than using Integrator<0, DimInt>::integrate_without_warning(double& value, double& error) const,
		 * which integrates over the unit hypercube. For maximal performance we thus recommend using that method instead
		 * and taking care of the integration boundaries by applying appropriate changes of variables in the integrand
		 * Integrator<0, DimInt>::Func.
		 */
		void integrate_without_warning (const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const;
		
		/**
		 * Assignment operator taking care of properly copying the integration Algorithm pointed to by Integrator<0, DimInt>::Alg
		 * from the Integrator<0, DimInt> \a otherIntegrator.
		 */
		Integrator& operator= (const Integrator& otherIntegrator);
		
		/**
		 * Destructor deleting the integration Algorithm pointed to by Integrator<0, DimInt>::Alg.
		 */
		~Integrator ();
		
	private:
		/**
		 * Function that shall be integrated.
		 */
		IntegrandWithoutFixedArguments<DimInt> Func;
		
		/**
		 * Pointer to the integration Algorithm.
		 */
		Algorithm* Alg;
		
		/**
		 * Lower integration boundaries.
		 * 
		 * These have to be \c mutable, as they will be modified by the methods
		 * Integrator<0, DimInt>::integrate(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const and
		 * Integrator<0, DimInt>::integrate_without_warning(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const.
		 */
		mutable Arguments<DimInt> LowerBounds;
		
		/**
		 * Upper integration boundaries.
		 * 
		 * These have to be \c mutable, as they will be modified by the methods
		 * Integrator<0, DimInt>::integrate(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const and
		 * Integrator<0, DimInt>::integrate_without_warning(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const.
		 */
		mutable Arguments<DimInt> UpperBounds;
		
		/**
		 * The \c string that will be written to the standard output by Integrator<0, DimInt>::error_handler to identify
		 * the specific Integrator<0, DimInt> object.
		 */
		std::string Identifier;
		
		/**
		 * Is called when the integration run performed by Integrator<0, DimInt>::integrate(double& value, double& error) const or
		 * Integrator<0, DimInt>::integrate(const Arguments<DimInt>& lowerBounds, const Arguments<DimInt>& uppperBounds, double& value, double& error) const
		 * fails. It writes Integrator<0, DimInt>::Identifier and all relevant information that can be extracted from the
		 * Algorithm::Result \a integral, i.e. the integral value, its estimated absolute and relative errors as well as
		 * possible further comments passed by the Algorithm pointed to by Integrator<0, DimInt>::Alg, to the standard output.
		 */
		void error_handler (const Algorithm::Result& integral) const;
		
		/**
		 * Wrapper for Integrator<0, DimInt>::Func that provides the form of the integrand expected by an integration Algorithm,
		 * i.e. an Algorithm::InternalIntegrand, if the integral shall be performed over the unit hypercube.
		 */
		double algorithm_internal_integrand_for_unit_hypercube (const double* dummyArgsFix, const double* argsInt) const;
		
		/**
		 * Wrapper for Integrator<0, DimInt>::Func that provides the form of the integrand expected by an integration Algorithm,
		 * i.e. an Algorithm::InternalIntegrand, if the integral shall be performed over the hypercube specified by \a LowerBounds
		 * and \a UpperBounds.
		 */
		double algorithm_internal_integrand_for_custom_hypercube (const double* dummyArgsFix, const double* argsInt) const;
	};
}

#include "Integrator.tpp"	// template implementations can not be compiled separately

#endif