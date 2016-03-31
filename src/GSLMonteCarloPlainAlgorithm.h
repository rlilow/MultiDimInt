#ifndef MULTIDIMINT_GSL_MONTE_CARLO_PLAIN_ALGORITHM_H
#define MULTIDIMINT_GSL_MONTE_CARLO_PLAIN_ALGORITHM_H

#include "GSLMonteCarloAlgorithm.h"

#include <string>

#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_rng.h>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a Monte Carlo integration scheme using the Plain algorithm of the GSL.
	 * 
	 * GSL documentation about <a href="https://www.gnu.org/software/gsl/manual/html_node/PLAIN-Monte-Carlo.html#PLAIN-Monte
	 * -Carlo">this algorithm</a>:
	 * 
	 * The plain Monte Carlo algorithm samples points randomly from the integration region to estimate the integral and
	 * its error.
	 */
	class GSLMonteCarloPlainAlgorithm : public GSLMonteCarloAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Plain algorithm of the GSL with absolute
		 * error limit \a absErr, relative error limit \a absRel, fixed number of integrand evaluations \a numEval, and
		 * GSL random number generator of type \a randomNumberGeneratorType.
		 *
		 * For further details on possible random number generator algorithms see the GSL documentation: <a href="https:
		 * //www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator
		 * -algorithms">Random number generator algorithms</a>.
		 */
		GSLMonteCarloPlainAlgorithm (double absErr, double relErr, std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType = gsl_rng_ranlxs2);
		
		GSLMonteCarloPlainAlgorithm* clone () const;
		
	protected:
		int gsl_mc_integration (std::size_t dimInt, gsl_monte_function& gslMonteIntegrand, double& value, double& error, std::string& furtherComment) const;
	};
}

#endif