#ifndef MULTIDIMINT_GSL_MONTE_CARLO_VEGAS_ALGORITHM_H
#define MULTIDIMINT_GSL_MONTE_CARLO_VEGAS_ALGORITHM_H

#include "GSLMonteCarloAlgorithm.hpp"

#include <string>

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a Monte Carlo integration scheme using the Vegas algorithm of the GSL.
	 * 
	 * GSL documentation about <a href="https://www.gnu.org/software/gsl/manual/html_node/VEGAS.html#VEGAS">this algorithm</a>:
	 * 
	 * The Vegas algorithm of Lepage is based on importance sampling. It samples points from the probability distribution
	 * described by the modulus of the integrand, so that the points are concentrated in the regions that make the largest
	 * contribution to the integral.
	 * 
	 * Author: Robert Lilow, ITA, ZAH, Heidelberg University (2016)
	 */
	class GSLMonteCarloVegasAlgorithm : public GSLMonteCarloAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Vegas algorithm of the GSL with absolute
		 * error limit \a absErr, relative error limit \a absRel, fixed number of integrand evaluations \a numEval, and
		 * GSL random number generator of type \a randomNumberGeneratorType. All further GSL Vegas specific parameters are
		 * set to the following standard values:
		 * 
		 *  - alpha = 1.5
		 *  - iterations = 5
		 *  - mode = 1 (GSL_VEGAS_MODE_IMPORTANCE)
		 * 
		 * For further details on these parameters or possible random number generator algorithms see the GSL documentation:
		 * <a href="https://www.gnu.org/software/gsl/manual/html_node/MISER.html#VEGAS">VEGAS</a>, <a href="https://www.gnu
		 * .org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator-algorithms">
		 * Random number generator algorithms</a>.
		 */
		GSLMonteCarloVegasAlgorithm (double absErr, double relErr, std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType = gsl_rng_ranlxs2);
		
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Vegas algorithm of the GSL with absolute
		 * error limit \a absErr, relative error limit \a absRel, fixed number of integrand evaluations \a numEval, and
		 * GSL random number generator of type \a randomNumberGeneratorType. All further arguments are GSL Vegas specific
		 * parameters.
		 * 
		 * For further details on these parameters or possible random number generator algorithms see the GSL documentation:
		 * <a href="https://www.gnu.org/software/gsl/manual/html_node/VEGAS.html#VEGAS">VEGAS</a>, <a href="https://www.gnu
		 * .org/software/gsl/manual/html_node/Random-number-generator-algorithms.html#Random-number-generator-algorithms">
		 * Random number generator algorithms</a>.
		 */
		GSLMonteCarloVegasAlgorithm (double absErr, double relErr, std::size_t numEval, const gsl_rng_type* randomNumberGeneratorType,
									 double alpha, std::size_t iterations, int mode);
		
		GSLMonteCarloVegasAlgorithm* clone () const;
		
	protected:
		int gsl_mc_integration (std::size_t dimInt, gsl_monte_function& gslMonteIntegrand, double& value, double& error, std::string& furtherComment) const;
		
	private:
		// GSL Vegas specific parameters
		double Alpha;
		int Iterations;
		int Mode;
	};
}

#endif