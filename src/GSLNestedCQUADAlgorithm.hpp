#ifndef MULTIDIMINT_GSL_NESTED_CQUAD_ALGORITHM_H
#define MULTIDIMINT_GSL_NESTED_CQUAD_ALGORITHM_H

#include "GSLNestedAlgorithm.hpp"

#include <gsl/gsl_math.h>

namespace MultiDimInt
{
	/**
	 * \brief Class implementing a nested integration scheme using the doubly-adaptive Clenshaw-Curtis-integration algorithm
	 * CQUAD of the GSL.
	 * 
	 * GSL documentation about <a href="https://www.gnu.org/software/gsl/manual/html_node/CQUAD-doubly_002dadaptive-integration
	 * .html#CQUAD-doubly_002dadaptive-integration">this algorithm</a>:
	 * 
	 * CQUAD is a new doubly-adaptive general-purpose quadrature routine which can handle most types of singularities,
	 * non-numerical function values such as Inf or NaN, as well as some divergent integrals. It generally requires more
	 * function evaluations than the integration routines in QUADPACK, yet fails less often for difficult integrands.
	 *
	 * The underlying algorithm uses a doubly-adaptive scheme in which Clenshaw-Curtis quadrature rules of increasing degree
	 * are used to compute the integral in each interval. The L_2-norm of the difference between the underlying interpolatory
	 * polynomials of two successive rules is used as an error estimate. The interval is subdivided if the difference
	 * between two successive rules is too large or a rule of maximum degree has been reached.
	 * 
	 * The CQUAD algorithm divides the integration region into subintervals, and in each iteration, the subinterval with
	 * the largest estimated error is processed. The algorithm uses Clenshaw-Curits quadrature rules of degree 4, 8, 16
	 * and 32 over 5, 9, 17 and 33 nodes respectively. Each interval is initialized with the lowest-degree rule. When an
	 * interval is processed, the next-higher degree rule is evaluated and an error estimate is computed based on the L_2-norm
	 * of the difference between the underlying interpolating polynomials of both rules. If the highest-degree rule has
	 * already been used, or the interpolatory polynomials differ significantly, the interval is bisected.
	 * 
	 * Author: Robert Lilow (2016)
	 */
	class GSLNestedCQUADAlgorithm : public GSLNestedAlgorithm
	{
	public:
		/**
		 * Constructor instantiating a nested GSL CQUAD integration algorithm with absolute error limit \a absErr, relative
		 * error limit \a relErr, and maximal number of stored intervals \a MaxInterval.
		 */
		GSLNestedCQUADAlgorithm (double absErr, double relErr, std::size_t maxInterval);
		
		GSLNestedCQUADAlgorithm* clone () const;
		
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