#ifndef MULTIDIMINT_ALGORITHM_H
#define MULTIDIMINT_ALGORITHM_H

#include <functional>
#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Abstract base class for integration algorithms.
	 * 
	 * Specific algorithms are implemented as classes derived from this base class.
	 * 
	 * Author: Robert Lilow, ITA, ZAH, Heidelberg University (2016)
	 */
	class Algorithm
	{
	public:
		/**
		 * Internal wrapper for the MultiDimInt::Integrand that shall be integrated which does not depend on the explicit 
		 * number of fixed arguments and integration variables, allowing integration algorithms to be independent of these
		 * as well. The meaning of the parameters stays the same, though. \a argsFix denotes the arguments which are kept
		 * fixed and \a argsInt contains the integration variables.
		 */
		using InternalIntegrand = std::function<double(const double* argsFix, const double* argsInt)>;
		
		/**
		 * Structure that contains all relevant results of an integration run. \a Failed is a \c bool that will be set to
		 * \c true if the integration failed, \a Value is the actual numerical value of the integral and \a Error the estimated
		 * absolute error. Furthermore, additional information can be passed via the \c string \a Comment.
		 * 
		 * These information are used by Integrator::integrate, Integrator:integrand_without_warning and Integrator::error_handler.
		 */
		struct Result
		{
			bool Failed;
			double Value;
			double Error;
			std::string Comment;
		};

		/**
		 * Performs the actual integration of the Algorithm::InternalIntegrand \a func using a specific algorithm and returns
		 * a Algorithm::Result \c struct containing all relevant results. \a argsFix are the fixed arguments that \a func
		 * depends on and \a dimInt is the number of its integration variables.
		 */
		virtual Result run (const InternalIntegrand& func, std::size_t dimInt, const double* argsFix) const = 0;
		
		/**
		 * Returns 'true' if this integration algorithm uses multiple cores in parallel and 'false' if it only runs
		 * on a single core.
		 * 
		 * This information can be used to decide if external paralleliziation is useful or unnecessary.
		 */
		virtual bool is_parallelized () const = 0;
		
		/**
		 * Dynamically creates a copy of this integration algorithm, using the copy-constructor, and returns a pointer
		 * to it.
		 * 
		 * This is needed to create a copy of some specific algorithm derived from this base class given only a pointer
		 * to a general Algorithm.
		 */
		virtual Algorithm* clone () const = 0;
		
		/**
		 * Default assignment operator.
		 */
		Algorithm& operator= (const Algorithm& otherAlgorithm) = default;
		
		/**
		 * Default destructor.
		 * 
		 * It is virtual to make sure that you can delete an instance of a specific algorithm derived from this base
		 * class through a pointer to a general Algorithm.
		 */
		virtual ~Algorithm () = default;
		
	protected:
		/**
		 * Constructor instantiating a general integration scheme with absolute error limit \a absErr and relative error
		 * limit \a relErr.
		 */
		Algorithm (double absErr, double relErr);
		
		/**
		 * Default copy-constructor.
		 */
		Algorithm (const Algorithm& otherAlgorithm) = default;
		
		/**
		 * Absolute error limit.
		 */
		double AbsErr;
		
		/**
		 * Relative error limit.
		 */
		double RelErr;
	};
}

#endif