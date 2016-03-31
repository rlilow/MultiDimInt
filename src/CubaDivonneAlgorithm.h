#ifndef MULTIDIMINT_CUBA_DIVONNE_ALGORITHM_H
#define MULTIDIMINT_CUBA_DIVONNE_ALGORITHM_H

#include "CubaAlgorithm.h"

#include <string>
#include <vector>

#include <cuba.h>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a Monte Carlo integration scheme using the Divonne algorithm of the Cuba library.
	 * 
	 * <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library documentation</a> about this algorithm:
	 * 
	 * Divonne uses stratified sampling for variance reduction, that is, it partitions the integration region such that
	 * all subregions have an approximately equal value of a quantity called the spread, defined as the product of the
	 * volume of the region and the difference between the maximum and the minimum of the integrand in this region. What
	 * sets Divonne apart from Suave is that the minimum and maximum of the integrand are sought using methods from numerical
	 * optimization. Particularly in high dimensions, the chance that one of the previously sampled points lies in or even
	 * close to the true extremum is fairly small.
	 * 
	 * On the other hand, the numerical minimization is beset with the usual pitfalls, i.e. starting from the lowest of a
	 * (relatively small) number of sampled points, Divonne will move directly into the local minimum closest to the starting
	 * point, which may or may not be close to the absolute minimum.
	 * 
	 * Divonne is a lot more complex than Suave and Vegas but also significantly faster for many integrands. 
	 */
	class CubaDivonneAlgorithm : public CubaAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Divonne algorithm of the Cuba library with
		 * absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations
		 * \a maxEval. All further Cuba Divonne specific parameters are set to the following values:
		 * 
		 *  - flags = 0
		 *  - mineval = 0
		 *  - statefile = \c NULL
		 *  - spin = \c NULL
		 *  - seed = 0
		 *  - key1 = 47
		 *  - key2 = 1
		 *  - key3 = 1
		 *  - maxpass = 5
		 *  - border = 0
		 *  - maxchisq = 10
		 *  - mindeviation = 0.25
		 *  - ngiven = 0
		 *  - xgiven = \c NULL
		 *  - nextra = 0
		 *  - peakfinder = \c NULL
		 * 
		 * For further details on these parameters see the <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library
		 * documentation</a>.
		 */
		CubaDivonneAlgorithm (double absErr, double relErr, int maxEval);
		
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using the Divonne algorithm of the Cuba library with
		 * absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations 
		 * \a maxEval. All further arguments are Cuba Divonne specific parameters.
		 * 
		 * For further details on these parameters see the <a href="http://arxiv.org/pdf/hep-ph/0404043.pdf">Cuba library
		 * documentation</a>.
		 */
		CubaDivonneAlgorithm (double absErr, double relErr, int maxEval,
							  int flags, int mineval, const std::string& statefile, void* spin,
							  int seed, int key1, int key2, int key3, int maxpass, double border, double maxchisq, double mindeviation, int ngiven, const std::vector<double>& xgiven, int nextra, peakfinder_t peakfinder);
		
		CubaDivonneAlgorithm* clone () const;
		
	protected:
		bool cuba_integration (int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const;
		
	private:
		// Cuba Divonne specific parameters
		int Seed;
		int Key1;
		int Key2;
		int Key3;
		int Maxpass;
		double Border;
		double Maxchisq;
		double Mindeviation;
		int Ngiven;
		std::vector<double> Xgiven;
		int Nextra;
		peakfinder_t Peakfinder;
	};
}

#endif