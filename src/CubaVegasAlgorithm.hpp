#ifndef MULTIDIMINT_CUBA_VEGAS_ALGORITHM_H
#define MULTIDIMINT_CUBA_VEGAS_ALGORITHM_H

#include "CubaAlgorithm.hpp"

#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a Monte Carlo integration scheme using the Vegas algorithm of the Cuba library.
	 * 
	 * <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library documentation</a> about this algorithm:
	 * 
	 * Vegas is a Monte Carlo algorithm that uses importance sampling as a variance-reduction technique. Vegas iteratively
	 * builds up a piecewise constant weight function, represented on a rectangular grid. Each iteration consists of a
	 * sampling step followed by a refinement of the grid.
	 * 
	 * Author: Robert Lilow (2016)
	 */
	class CubaVegasAlgorithm : public CubaAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Vegas algorithm of the Cuba library with
		 * absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations
		 * \a maxEval. All further Cuba Vegas specific parameters are set to the following values:
		 * 
		 *  - flags = 0
		 *  - mineval = 0
		 *  - statefile = ""
		 *  - spin = \c NULL
		 *  - seed = 0
		 *  - nstart = 1000
		 *  - nincrease = 500
		 *  - nbatch = 1000
		 *  - gridno = 0
		 * 
		 * For further details on these parameters see the <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library
		 * documentation</a>.
		 */
		CubaVegasAlgorithm (double absErr, double relErr, int maxEval);
		
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Vegas algorithm of the Cuba library with
		 * absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations 
		 * \a maxEval. All further arguments are Cuba Vegas specific parameters.
		 * For further details on these parameters see the <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library
		 * documentation</a>.
		 */
		CubaVegasAlgorithm (double absErr, double relErr, int maxEval,
							int flags, int mineval, const std::string& statefile, void* spin,
							int seed, int nstart, int nincrease, int nbatch, int gridno);
		
		CubaVegasAlgorithm* clone () const;
		
	protected:
		bool cuba_integration (int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const;
		
	private:
		// Cuba Vegas specific parameters
		int Seed;
		int Nstart;
		int Nincrease;
		int Nbatch;
		int Gridno;
	};
}

#endif