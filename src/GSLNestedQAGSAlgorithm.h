#ifndef MULTIDIMINT_GSL_NESTED_QAGS_ALGORITHM_H
#define MULTIDIMINT_GSL_NESTED_QAGS_ALGORITHM_H

#include "GSLNestedAlgorithm.h"

#include <gsl/gsl_math.h>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a nested integration scheme using the adaptive Gauss-Kronrod-integration algorithm with
	 * singularities QAGSS of the GSL.
	 * 
	 * GSL documentation about <a href="https://www.gnu.org/software/gsl/manual/html_node/QAGSS-adaptive-integration-with
	 * -singularities.html#QAGSS-adaptive-integration-with-singularities">this algorithm</a>:
	 * 
	 * The presence of an integrable singularity in the integration region causes an adaptive routine to concentrate new
	 * subintervals around the singularity. As the subintervals decrease in size the successive approximations to the integral
	 * converge in a limiting fashion. This approach to the limit can be accelerated using an extrapolation procedure. The
	 * QAGSS algorithm combines adaptive bisection with the Wynn epsilon-algorithm to speed up the integration of many types
	 * of integrable singularities. 
	 */
	class GSLNestedQAGSAlgorithm : public GSLNestedAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a nested GSL QAGS integration algorithm with absolute error limit \a absErr, relative
		 * error limit \a relErr, and maximal number of intervals \a MaxInterval.
		 */
		GSLNestedQAGSAlgorithm (double absErr, double relErr, std::size_t maxInterval);
		
		GSLNestedQAGSAlgorithm* clone () const;
		
	protected:
		int gsl_integration (const gsl_function& recursiveIntegrand, double& value, double& error) const;
		
	private:
		/**
		 * Maximal number of intervals.
		 */
		std::size_t MaxInterval;
	};
}

#endif