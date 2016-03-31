#ifndef MULTIDIMINT_GSL_MONTE_CARLO_MISER_ALGORITHM_H
#define MULTIDIMINT_GSL_MONTE_CARLO_MISER_ALGORITHM_H

#include "GSLMonteCarloAlgorithm.h"

#include <string>

#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_rng.h>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a Monte Carlo integration scheme using the Miser algorithm of the GSL.
	 * 
	 * GSL documentation about <a href="https://www.gnu.org/software/gsl/manual/html_node/MISER.html#MISER">this algorithm</a>:
	 * 
	 * The Miser algorithm of Press and Farrar is based on recursive stratified sampling. This technique aims to reduce
	 * the overall integration error by concentrating integration points in the regions of highest variance. 
	 */
	class GSLMonteCarloMiserAlgorithm : public GSLMonteCarloAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Miser algorithm of the GSL with absolute
		 * error limit \a absErr, relative error limit \a absRel, fixed number of integrand evaluations \a numEval, and
		 * GSL random number generator of type \a randomNumberGeneratorType. All further GSL Miser specific parameters are
		 * set to the following standard values:
		 * 
		 *  - estimate_frac = 0.1
		 *  - min_calls_per_dim = 16
		 *  - min_calls_per_dim_per_bisection = 32 * min_calls_per_dim = 512
		 *  - alpha = 2
		 *  - dither = 0
		 * 
		 * For further details on these parameters or possible random number generator algorithms see the GSL documentation:
		 * <a href="https://www.gnu.org/software/gsl/manual/html_node/MISER.html#MISER">MISER</a>, <a href="https://www.gnu
		 * .org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator-algorithms">
		 * Random number generator algorithms</a>.
		 */
		GSLMonteCarloMiserAlgorithm (double absErr, double relErr, std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType = gsl_rng_ranlxs2);
		
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Miser algorithm of the GSL with absolute
		 * error limit \a absErr, relative error limit \a absRel, fixed number of integrand evaluations \a numEval, and
		 * GSL random number generator of type \a randomNumberGeneratorType. All further arguments are GSL Miser specific
		 * parameters.
		 * 
		 * For further details on these parameters or possible random number generator algorithms see the GSL documentation:
		 * <a href="https://www.gnu.org/software/gsl/manual/html_node/MISER.html#MISER">MISER</a>, <a href="https://www.gnu
		 * .org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator-algorithms">
		 * Random number generator algorithms</a>.
		 */
		GSLMonteCarloMiserAlgorithm (double absErr, double relErr, std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType,
									 double estimate_frac, std::size_t min_calls_per_dim, std::size_t min_calls_per_dim_per_bisection, double alpha, double dither);
		
		GSLMonteCarloMiserAlgorithm* clone () const;
		
	protected:
		int gsl_mc_integration (std::size_t dimInt, gsl_monte_function& gslMonteIntegrand, double& value, double& error, std::string& furtherComment) const;
		
	private:
		// GSL Miser specific parameters
		double Estimate_frac;
		int Min_calls_per_dim;
		int Min_calls_per_dim_per_bisection;
		double Alpha;
		double Dither;
	};
}

#endif