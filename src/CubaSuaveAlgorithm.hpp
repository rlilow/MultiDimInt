#ifndef MULTIDIMINT_CUBA_SUAVE_ALGORITHM_H
#define MULTIDIMINT_CUBA_SUAVE_ALGORITHM_H

#include "CubaAlgorithm.hpp"

#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a Monte Carlo integration scheme using the Suave algorithm of the Cuba library.
	 * 
	 * <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library documentation</a> about this algorithm:
	 * 
	 * Suave (short for subregion-adaptive vegas) uses Vegas-like importance sampling combined with a globally adaptive
	 * subdivision strategy: Until the requested accuracy is reached, the region with the largest error at the time is bisected
	 * in the dimension in which the fluctuations of the integrand are reduced most. The number of new samples in each half
	 * is prorated for the fluctuation in that half.
	 * 
	 * A similar method, known as recursive stratified sampling, is implemented in Miser. Miser always samples a fixed number
	 * of points, however, which is somewhat undesirable since it does not stop once the prescribed accuracy is reached.
	 * 
	 * Suave first samples the integration region in a Vegas-like step, i.e. using importance sampling with a separable
	 * weight function. It then slices the integration region in two, as Miser would do. Suave does not immediately recurse
	 * on those subregions, however, but maintains a list of all subregions and selects the region with the largest absolute
	 * error for the next cycle of sampling and subdivision. That is, Suave uses global error estimation and terminates
	 * when the requested relative or absolute accuracy is attained.
	 * 
	 * Author: Robert Lilow (2016)
	 */
	class CubaSuaveAlgorithm : public CubaAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Suave algorithm of the Cuba library with
		 * absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations
		 * \a maxEval. All further Cuba Suave specific parameters are set to the following values:
		 * 
		 *  - flags = 0
		 *  - mineval = 0
		 *  - statefile = ""
		 *  - spin = \c NULL
		 *  - seed = 0
		 *  - nnew = 1000
		 *  - nmin = 2
		 *  - flatness = 50
		 * 
		 * For further details on these parameters see the <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library
		 * documentation</a>.
		 */
		CubaSuaveAlgorithm (double absErr, double relErr, int maxEval);
		
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Suave algorithm of the Cuba library with
		 * absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations 
		 * \a maxEval. All further arguments are Cuba Suave specific parameters.
		 * 
		 * For further details on these parameters see the <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library
		 * documentation</a>.
		 */
		CubaSuaveAlgorithm (double absErr, double relErr, int maxEval,
							int flags, int mineval, const std::string& statefile, void* spin,
							int seed, int nnew, int nmin, double flatness);
		
		CubaSuaveAlgorithm* clone () const;
		
	protected:
		bool cuba_integration (int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const;
		
	private:
		// Cuba Suave specific parameters
		int Seed;
		int Nnew;
		int Nmin;
		double Flatness;
	};
}

#endif