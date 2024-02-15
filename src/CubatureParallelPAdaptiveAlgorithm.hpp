#ifndef MULTIDIMINT_CUBATURE_PARALLEL_PADAPTIVE_ALGORITHM_H
#define MULTIDIMINT_CUBATURE_PARALLEL_PADAPTIVE_ALGORITHM_H

#include "CubatureParallelAlgorithm.hpp"

#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a parallelised integration algorithm using the vectorised p-adaptive cubature algorithm
	 * of the <a href="https://github.com/stevengj/cubature">Cubature libarary</a> by Steven G. Johnson.
	 *
	 * Author: Robert Lilow (2024)
	 */
	class CubatureParallelPAdaptiveAlgorithm : public CubatureParallelAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a parallel integration scheme using the vectorised p-adaptive cubature algorithm of
		 * the Cubature library with absolute error limit \a absErr, relative error limit \a absRel and maximal number
		 * of integrand evaluations \a maxEval.
		 */
		CubatureParallelPAdaptiveAlgorithm(double absErr, double relErr, std::size_t maxEval);

		CubatureParallelPAdaptiveAlgorithm *clone() const;

	protected:
		bool cubature_integration(unsigned dimInt, CubatureData &cubatureData, double &value, double &error, std::string &furtherComment) const;
	};
}

#endif
