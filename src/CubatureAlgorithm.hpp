#ifndef MULTIDIMINT_CUBATURE_ALGORITHM_H
#define MULTIDIMINT_CUBATURE_ALGORITHM_H

#include "Algorithm.hpp"

#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Abstract base class for integration algorithms using the <a
	 * href="https://github.com/stevengj/cubature">Cubature libarary</a> by Steven G. Johnson.
	 *
	 * The Cubature library provides two different cubature integration schemes, both in a standard and a vectorised
	 * form (hcubature, pcubature, hcubature_v, pcubature_v).
	 *
	 * Author: Robert Lilow (2024)
	 */
	class CubatureAlgorithm : public Algorithm
	{
	public:
		Algorithm::Result run(const InternalIntegrand &func, std::size_t dimInt, const double *argsFix) const;

		virtual bool is_parallelized() const = 0;

		virtual CubatureAlgorithm *clone() const = 0;

	protected:
		/**
		 * Constructor instantiating an integration scheme using some algorithm of the Cubature library with absolute
		 * error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations \a maxEval.
		 */
		CubatureAlgorithm(double absErr, double relErr, std::size_t maxEval);

		/**
		 * Structure gathering all information needed by CubatureAlgorithm::cubature_integration and
		 * CubatureAlgorithm::cubature_integrand. \a Func is the Algorithm::InternalIntegrand to be integrated for fixed
		 * arguments \a ArgsFix.
		 */
		struct CubatureData
		{
			CubatureData(const InternalIntegrand &func, const double *argsFix) : Func(func),
																				 ArgsFix(argsFix){};

			const InternalIntegrand &Func;
			const double *ArgsFix;
		};

		/**
		 * Performs the actual integration of CubatureAlgorithm::cubature_integrand by calling the appropriate function
		 * of the Cubature library. It takes the number of integration variables \a dimInt, a reference to a
		 * CubatureAlgorithm::CubatureData \c struct \a cubatureData provided by CubatureAlgorithm::run, writes the
		 * numerical value of the integral into \a value, the estimated error into \a error and an optional further
		 * comment into \a furtherComment. It returns \c true if the integration succeeded and \c false otherwise.
		 */
		virtual bool cubature_integration(unsigned dimInt, CubatureData &cubatureData, double &value, double &error, std::string &furtherComment) const = 0;

		/**
		 * Maximal number of integrand evaluations.
		 */
		int MaxEval;
	};
}

#endif
