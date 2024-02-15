#ifndef MULTIDIMINT_CUBATURE_SERIAL_HADAPTIVE_ALGORITHM_H
#define MULTIDIMINT_CUBATURE_SERIAL_HADAPTIVE_ALGORITHM_H

#include "CubatureSerialAlgorithm.hpp"

#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a serialised integration algorithm using the non-vectorised h-adaptive cubature
	 * algorithm of the <a href="https://github.com/stevengj/cubature">Cubature libarary</a> by Steven G. Johnson.
	 *
	 * Author: Robert Lilow (2024)
	 */
	class CubatureSerialHAdaptiveAlgorithm : public CubatureSerialAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a serial integration scheme using the non-vectorised h-adaptive cubature algorithm
		 * of the Cubature library with absolute error limit \a absErr, relative error limit \a absRel and maximal
		 * number of integrand evaluations \a maxEval.
		 */
		CubatureSerialHAdaptiveAlgorithm(double absErr, double relErr, std::size_t maxEval);

		CubatureSerialHAdaptiveAlgorithm *clone() const;

	protected:
		bool cubature_integration(unsigned dimInt, CubatureData &cubatureData, double &value, double &error, std::string &furtherComment) const;
	};
}

#endif
