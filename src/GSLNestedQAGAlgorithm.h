#ifndef MULTIDIMINT_GSL_NESTED_QAG_ALGORITHM_H
#define MULTIDIMINT_GSL_NESTED_QAG_ALGORITHM_H

#include "GSLNestedAlgorithm.h"

#include <gsl/gsl_math.h>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a nested integration scheme using the adaptive Gauss-Kronrod-integration algorithm QAG
	 * of the GSL.
	 * 
	 * GSL documentation about <a href="https://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html#QAG
	 * -adaptive-integration">this algorithm</a>:
	 * 
	 * The QAG algorithm is a simple adaptive integration procedure. The integration region is divided into subintervals,
	 * and on each iteration the subinterval with the largest estimated error is bisected. This reduces the overall error
	 * rapidly, as the subintervals become concentrated around local difficulties in the integrand.
	 * 
	 * Author: Robert Lilow, ITA, ZAH, Heidelberg University (2016)
	 */
	class GSLNestedQAGAlgorithm : public GSLNestedAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a nested GSL QAG integration algorithm with absolute error limit \a absErr, relative
		 * error limit \a relErr, maximal number of intervals \a MaxInterval, and Gauss-Kronrod rule \a key. The \a key
		 * value has to lie in between 1 and 6, corresponding to a 15, 21, 31, 41, 51 and 61 point rule.
		 */
		GSLNestedQAGAlgorithm (double absErr, double relErr, std::size_t maxInterval, int key = 1);
		
		GSLNestedQAGAlgorithm* clone () const;
		
	protected:
		int gsl_integration (const gsl_function& recursiveIntegrand, double& value, double& error) const;
		
	private:
		/**
		 * Maximal number of intervals.
		 */
		std::size_t MaxInterval;
		
		/**
		 * Gauss-Kronrod rule.
		 */
		int Key;
	};
}

#endif