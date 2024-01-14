#ifndef MULTIDIMINT_CUBA_ALGORITHM_H
#define MULTIDIMINT_CUBA_ALGORITHM_H

#include "Algorithm.h"

#include <string>

namespace MultiDimInt
{
	/**
	 * \brief Abstract base class for integration algorithms of the Cuba libarary.
	 * 
	 * The Cuba library provides three different Monte Carlo integration schemes (Vegas, Suave, Divonne) and one cubature
	 * scheme (Cuhre). This class acts as a wrapper for these algorithms.
	 * 
	 * Author: Robert Lilow, ITA, ZAH, Heidelberg University (2016)
	 */
	class CubaAlgorithm : public Algorithm
	{
	public:
		Algorithm::Result run (const InternalIntegrand& func, std::size_t dimInt, const double* argsFix) const;
		
		bool is_parallelized () const;
		
		virtual CubaAlgorithm* clone () const = 0;
		
	protected:
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using some algorithm of the Cuba library with absolute
		 * error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations \a maxEval.
		 * The common Cuba parameters are set to the following values:
		 * 
		 *  - flags = 0
		 *  - mineval = 0
		 *  - statefile = ""
		 *  - spin = \c NULL
		 *  - key = 0
		 * 
		 * For further details on these parameters or possible specific Cuba algorithms see the <a href="http://arxiv.org
		 * /pdf/hep-ph/0404043.pdf">Cuba library documentation</a>.
		 */
		CubaAlgorithm (double absErr, double relErr, int maxEval);
		
		/**
		 * Constructor instantiating a Monte Carlo integration scheme using some algorithm of the Cuba library with absolute
		 * error limit \a absErr, relative error limit \a absRel and maximal number of integrand evaluations \a maxEval.
		 * All further arguments are the common Cuba parameters.
		 *
		 * For further details on these parameters or possible specific Cuba algorithms see the <a href="http://arxiv.org
		 * /pdf/hep-ph/0404043.pdf">Cuba library documentation</a>.
		 */
		CubaAlgorithm (double absErr, double relErr, int maxEval,
					   int flags, int mineval, const std::string& statefile, void* spin);
		
		/**
		 * Structure gathering all information needed by CubaAlgorithm::cuba_integration and CubaAlgorithm::cuba_integrand.
		 * \a Func is the Algorithm::InternalIntegrand to be integrated for fixed arguments \a ArgsFix.
		 */
		struct CubaData
		{
			CubaData (const InternalIntegrand& func, const double* argsFix) :
				Func(func),
				ArgsFix(argsFix)
			{};
			
			const InternalIntegrand& Func;
			const double* ArgsFix;
		};
		
		/**
		 * Performs the actual integration of CubaAlgorithm::cuba_integrand by calling the appropriate function of the Cuba
		 * library. It takes the number of integration variables \a dimInt, a reference to a CubaAlgorithm::CubaData \c struct
		 * \a cubaData provided by CubaAlgorithm::run, writes the numerical value of the integral into \a value, the estimated
		 * error into \a error, the probability that this error estimate is wrong into \a prop and an optional further comment
		 * into \a furtherComment. It returns \c true if the integration succeeded and \c false otherwise.
		 */
		virtual bool cuba_integration (int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const = 0;
		
		/**
		 * Wrapper for the function to be integrated that provides the form of the integrand expected by Cuba integration
		 * routines. It has to be \c static, as the Cuba integration routines only accept non-member functions. Therefore,
		 * all information needed by the integrand has to be provided via the \c void pointer \a cubaData, since static
		 * methods cannot access the non-public members of their class. \a dimInt is the number of integration variables
		 * gathered in \a argsInt, \a ncomp is the number of components of the integrand (always 1 in MultiDimInt::Function),
		 * and the value of the integral is written into \a value.
		 */
		static int cuba_integrand (const int* dimInt, const double* argsInt, const int* ncomp, double* value, void* cubaData);
		
		/**
		 * Maximal number of integrand evaluations.
		 */
		int MaxEval;
		
		// Cuba specific parameters
		int Flags;
		int Mineval;
		std::string Statefile;
		void* Spin;
		int Key;
	};
}

#endif
