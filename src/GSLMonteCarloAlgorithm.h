#ifndef MULTIDIMINT_GSL_MONTE_CARLO_ALGORITHM_H
#define MULTIDIMINT_GSL_MONTE_CARLO_ALGORITHM_H

#include "Algorithm.h"

#include <string>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

namespace MultiDimInt
{
	/**
	 * \brief Abstract base class for Monte Carlo integration algorithms of the GSL.
	 * 
	 * The GSL provides three different Monte Carlo integration schemes (Plain, Miser, Vegas). This class acts as a wrapper
	 * for these algorithms.
	 * 
	 * Author: Robert Lilow, ITA, ZAH, Heidelberg University (2016)
	 */
	class GSLMonteCarloAlgorithm : public Algorithm
	{
	public:
		Algorithm::Result run (const InternalIntegrand& func, std::size_t dimInt, const double* argsFix) const;
		
		bool is_parallelized () const;
		
		virtual GSLMonteCarloAlgorithm* clone () const = 0;
		
		/**
		 * Assignment operator taking care of properly copying the GSL random number generator GSLMonteCarloAlgorithm::RandomNumberGeneratorType
		 * from the Algorithm \a otherGSLMonteCarloAlgorithm.
		 */
		GSLMonteCarloAlgorithm& operator= (const GSLMonteCarloAlgorithm& otherGSLMonteCarloAlgorithm);
		
		/**
		 * Destructor freeing the GSL random number generator GSLMonteCarloAlgorithm::RandomNumberGenerator.
		 * 
		 * It is virtual to make sure that you can delete an instance of a specific algorithm derived from this base
		 * class through a pointer to a general Algorithm.
		 */
		virtual ~GSLMonteCarloAlgorithm ();
		
	protected:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using some algorithm of the GSL with absolute
		 * error limit \a absErr, relative error limit \a absRel, fixed number of integrand evaluations \a numEval, and
		 * GSL random number generator of type \a randomNumberGeneratorType.
		 *
		 * For further details on possible specific Monte Carlo algorithms or possible random number generator algorithms
		 * see the GSL documentation: <a href="https://www.gnu.org/software/gsl/manual/html_node/Monte-Carlo-Integration
		 * .html#Monte-Carlo-Integration">Monte Carlo Integration</a>, <a href="https://www.gnu.org/software/gsl/manual
		 * /html_node/Random-number-generator-algorithms.html#Random-number* -generator-algorithms">Random number generator
		 * algorithms</a>.
		 */
		GSLMonteCarloAlgorithm (double absErr, double relErr, std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType = gsl_rng_ranlxs2);
		
		/**
		 * Copy-constructor taking care of properly copying the GSL random number generator GSLMonteCarloAlgorithm::RandomNumberGeneratorType
		 * from the Algorithm \a otherGSLMonteCarloAlgorithm.
		 */
		GSLMonteCarloAlgorithm (const GSLMonteCarloAlgorithm& otherGSLMonteCarloAlgorithm);
		
		/**
		 * Structure gathering all information needed by GSLMonteCarloAlgorithm::gsl_mc_integrand. \a Func is the
		 * Algorithm::InternalIntegrand to be integrated for fixed arguments \a ArgsFix.
		 */
		struct GSLMonteCarloData
		{
			GSLMonteCarloData (const InternalIntegrand& func, const double* argsFix) :
				Func(func),
				ArgsFix(argsFix)
			{};
			
			const InternalIntegrand& Func;
			const double* ArgsFix;
		};
		
		/**
		 * Performs the actual integration of GSLMonteCarloAlgorithm::gsl_mc_integrand by calling the appropriate function
		 * of the GSL. \a dimInt is the number of integration variables and \a gslMonteIntegrand a pointer to the \c gsl_monte_function
		 * expected by the integration routine. It writes the numerical value of the integral into \a value, the estimated
		 * error into \a error, the probability that this error estimate is wrong into \a prop and an optional further comment
		 * into \a furtherComment. Furthermore, it returns the GSL error code.
		 */
		virtual int gsl_mc_integration (std::size_t dimInt, gsl_monte_function& gslMonteIntegrand, double& value, double& error, std::string& furtherComment) const = 0;
		
		/**
		 * Wrapper for the function to be integrated that provides the form of the integrand expected by GSL Monte Carlo
		 * integration routines. It has to be \c static, as the GSL Monte Carlo integration routines only accept non-member
		 * functions. Therefore, all information needed by the integrand has to be provided via the \c void pointer
		 * \a gslMonteCarloData, since static methods cannot access the non-public members of their class. \a dimInt is the
		 * number of integration variables gathered in \a argsInt. The value of the integral is returned.
		 */
		static double gsl_mc_integrand (double* argsInt, std::size_t dimInt, void* gslMonteCarloData);
		
		/**
		 * Fixed number of integrand evaluations.
		 */
		std::size_t NumEval;
		
		/**
		 * Type of the GSL random number generator used for the sampling.
		 */
		const gsl_rng_type* RandomNumberGeneratorType;
		
		/**
		 * GSL random number generator used for the sampling.
		 */
		gsl_rng* RandomNumberGenerator;
	};
}

#endif