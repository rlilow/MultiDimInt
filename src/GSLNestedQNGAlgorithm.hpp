#ifndef MULTIDIMINT_GSL_NESTED_QNG_ALGORITHM_H
#define MULTIDIMINT_GSL_NESTED_QNG_ALGORITHM_H

#include "GSLNestedAlgorithm.hpp"

#include <gsl/gsl_math.h>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a nested integration scheme using the non-adaptive Gauss-Kronrod-integration algorithm QNG
	 * of the GSL.
	 * 
	 * GSL documentation about <a href="https://www.gnu.org/software/gsl/manual/html_node/QNG-non_002dadaptive-Gauss_002dKronrod
	 * -integration.html#QNG-non_002dadaptive-Gauss_002dKronrod-integration">this algorithm</a>:
	 * 
	 * The QNG algorithm is a non-adaptive procedure which uses fixed Gauss-Kronrod-Patterson abscissae to sample the integrand
	 * at a maximum of 87 points. It is provided for fast integration of smooth functions.
	 * 
	 * Author: Robert Lilow (2016)
	 */
	class GSLNestedQNGAlgorithm : public GSLNestedAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a nested GSL QNG integration algorithm with absolute error limit \a absErr and relative
		 * error limit \a relErr. It does not take a further argument specifying the maximal number of internal integration
		 * steps, as the GSL QNG algorithm uses fixed numbers of sampling points.
		 */
		GSLNestedQNGAlgorithm (double absErr, double relErr);
		
		GSLNestedQNGAlgorithm* clone () const;
		
	protected:
		int gsl_integration (const gsl_function& recursiveIntegrand, double& value, double& error) const;
	};
}

#endif