#ifndef MULTIDIMINT_CUBA_CUHRE_ALGORITHM_H
#define MULTIDIMINT_CUBA_CUHRE_ALGORITHM_H

#include "CubaAlgorithm.h"

#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a Monte Carlo integration scheme using the Cuhre algorithm of the Cuba library.
	 * 
	 * <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library documentation</a> about this algorithm:
	 * 
	 * Cuhre is a deterministic algorithm which uses one of several cubature rules of polynomial degree in a globally adaptive
	 * subdivision scheme. The subdivision algorithm is similar to Suave's.
	 * 
	 * In moderate dimensions Cuhre is very competitive, particularly if the integrand is well approximated by polynomials.
	 * As the dimension increases, the number of points sampled by the cubature rules rises considerably, however, and by
	 * the same token the usefulness declines.
	 */
	class CubaCuhreAlgorithm : public CubaAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Cuhre algorithm of the Cuba library with
		 * absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations
		 * \a maxEval. The Cuba Cuhre specific parameters are set to the following values:
		 * 
		 *  - flags = 0
		 *  - mineval = 0
		 *  - statefile = \c NULL
		 *  - spin = \c NULL
		 *  - key = 0
		 * 
		 * For further details on these parameters see the <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library
		 * documentation</a>.
		 */
		CubaCuhreAlgorithm (double absErr, double relErr, int maxEval);
		
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Cuhre algorithm of the Cuba library with
		 * absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations 
		 * \a maxEval. All further arguments are Cuba Cuhre specific parameters.
		 *
		 * For further details on these parameters see the <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library
		 * documentation</a>.
		 */
		CubaCuhreAlgorithm (double absErr, double relErr, int maxEval,
							int flags, int mineval, const std::string& statefile, void* spin,
							int key);
		
		CubaCuhreAlgorithm* clone () const;
		
	protected:
		bool cuba_integration (int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const;
		
	private:
		// Cuba Cuhre specific parameter
		int Key;
	};
}

#endif