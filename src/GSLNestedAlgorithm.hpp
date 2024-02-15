#ifndef MULTIDIMINT_GSL_NESTED_ALGORITHM_H
#define MULTIDIMINT_GSL_NESTED_ALGORITHM_H

#include "Algorithm.hpp"

#include <vector>

#include <gsl/gsl_math.h>

namespace MultiDimInt
{
	/**
	 * \brief Abstract base class for nested integration schemes using integration algorithms of the GSL.
	 * 
	 * The nesting of one-dimensional GSL integrations is implemented via a recursion procedure. Specific algorithms are
	 * implemented as classes derived from this base class. They are expected to use one of the algorithms listed in the
	 * GSL documentation under <a href="https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html#Numerical
	 * -Integration">Numerical Integration</a>.
	 * 
	 * Note that the error handling is only performed for the outermost integration.
	 * 
	 * Author: Robert Lilow (2016)
	 */
	class GSLNestedAlgorithm : public Algorithm
	{
	public:
		Algorithm::Result run (const InternalIntegrand& func, std::size_t dimInt, const double* argsFix) const;
		
		bool is_parallelized () const;
		
		virtual GSLNestedAlgorithm* clone () const = 0;
		
	protected:
		/**
		 * Constructor instantiating a general nested GSL integration algorithm with absolute error limit \a absErr and relative
		 * error limit \a relErr.
		 */
		GSLNestedAlgorithm (double absErr, double relErr);
		
		/**
		 * Implements the nested integration by recursively calling itself. Each recursion step treats a single one-dimensional
		 * integration over the argument \a currentIntArg. This function is of the form expected for the integrand in a
		 * gsl_function \c struct. It has to be static, as the GSL integration routines only accept non-member functions.
		 * Therefore, all information needed by each recursion step has to be provided via the \c void pointer \a recursionData,
		 * since static methods cannot access the non-public members of their class. It expects a GSLNestedAlgorithm::GSLNestedData
		 * pointer.
		 */
		static double internal_recursion (double currentIntArg, void* gslNestedData);
		
		/**
		 * Performs the actual one-dimensional GSL integration on the \c gsl_function \a recursiveIntegrand that will be
		 * provided by GSLNestedAlgorithm::internal_recursion. It writes the numerical value of the integral into \a value
		 * and the estimated error into \a error, and returns the GSL error code.
		 */
		virtual int gsl_integration (const gsl_function& recursiveIntegrand, double& value, double& error) const = 0;
		
		/**
		 * Structure gathering all information needed in each recursion step of GSLNestedAlgorithm::internal_recursion.
		 * \a NestingCounter keeps track of the current recursion depth, \a Integrand is the Algorithm::AlgorithmFunction
		 * to be integrated for fixed arguments \a ArgsFix, \a ArgsInt contains the current values of the \a DimInt integration
		 * variables, \a ThisGSLNestedAlgorithm is a pointer to the GSLNestedAlgorithm object itself, and \a Integral gathers
		 * all relevant information of the integration to be returned by GSLNestedAlgorithm::run.
		 */
		struct GSLNestedData
		{
			GSLNestedData (std::size_t nestingCounter, const InternalIntegrand& func, const double* argsFix, double* argsInt, std::size_t dimInt, const GSLNestedAlgorithm* thisGSLNestedAlgorithm) :
				NestingCounter(nestingCounter),
				Func(func),
				ArgsFix(argsFix),
				ArgsInt(argsInt),
				DimInt(dimInt),
				ThisGSLNestedAlgorithm(thisGSLNestedAlgorithm)
			{};
			
			std::size_t NestingCounter;
			const InternalIntegrand& Func;
			const double* ArgsFix;
			double* ArgsInt;
			const std::size_t DimInt;
			const GSLNestedAlgorithm* ThisGSLNestedAlgorithm;
		};
	};
}

#endif