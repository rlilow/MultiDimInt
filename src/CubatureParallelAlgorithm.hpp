#ifndef MULTIDIMINT_CUBATURE_PARALLEL_ALGORITHM_H
#define MULTIDIMINT_CUBATURE_PARALLEL_ALGORITHM_H

#include "CubatureAlgorithm.hpp"

#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Abstract base class for parallelised integration algorithms using the <a
	 * href="https://github.com/stevengj/cubature">Cubature libarary</a> by Steven G. Johnson.
	 *
	 * This class uses the two vectorised cubature integration schemes (hcubature_v, pcubature_v) to parallelise calls
	 * to the integrand.
	 *
	 * Author: Robert Lilow (2024)
	 */
	class CubatureParallelAlgorithm : public CubatureAlgorithm
	{
	public:
		bool is_parallelized() const;

		virtual CubatureParallelAlgorithm *clone() const = 0;

	protected:
		/**
		 * Constructor instantiating a parallel integration scheme using some vectorised algorithm of the Cubature
		 * library with absolute error limit \a absErr, relative error limit \a absRel and maximal number of integrand
		 * evaluations \a maxEval.
		 */
		CubatureParallelAlgorithm(double absErr, double relErr, std::size_t maxEval);

		virtual bool cubature_integration(unsigned dimInt, CubatureData &cubatureData, double &value, double &error, std::string &furtherComment) const = 0;

		/**
		 * Wrapper for the function to be integrated that provides the form of the integrand expected by vectorized
		 * Cubature integration routines. It has to be \c static, as the Cubature integration routines only accept
		 * non-member functions. Therefore, all information needed by the integrand has to be provided via the \c void
		 * pointer \a cubatureData, since static methods cannot access the non-public members of their class. \a dimInt
		 * is the number of integration variables gathered in \a argsInt, \a dimFunc is the number of components of the
		 * integrand (always 1 in MultiDimInt::Function), \a numPoints is the number of points to be evaluated in
		 * parallel, and the value of the integral is written into \a value.
		 */
		static int cubature_integrand(unsigned dimInt, std::size_t numPoints, const double *argsInt, void *cubatureData, unsigned dimFunc, double *value);
	};
}

#endif
