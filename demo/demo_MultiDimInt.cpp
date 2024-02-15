#include "../MultiDimInt.hpp"

#include <cmath>
#include <iostream>

/**
 * MultiDimInt demo:
 * 
 * Integral of the function f(x;y) = f(x0,x1;y0,y1,y2) = x0*y0 * (x1-y1)^4 * exp(-y2^2) over its arguments y0, y1 and y2
 * for fixed arguments x0=2 and x1=1.
 * 
 * MultiDimInt allows to perform this integral for different finite as well as infinite integration boundaries, and also
 * provides many different integration algorithms. Additionally, it is able to handle various ways of implementing this
 * function.
 * 
 * To get an idea of all these different possibilities, 5 exemplary cases are demonstrated:
 * 
 *  1) boundaries:     0<y0<1, 0<y1<1, 0<y2<1
 *     algorithm:      nested GSL CQUAD (nested quadrature)
 *     implementation: global function that expects x and y as arguments
 * 
 *  2) boundaries:     0<y0<3, -1<y1<1, -2<y2<3
 *     algorithm:      GSL Vegas (Monte Carlo)
 *     implementation: function object that expects x and y as arguments
 * 
 *  3) boundaries:     0<y0<3, -1<y1<1, -2<y2<infinity
 *     algorithm:      Cuba Divonne (Monte Carlo)
 *     implementation: function object that expects only y as argument and sets x via the constructor
 * 
 *  4) boundaries:     0<y0<3, -1<y1<1, -infinity<y2<infinity
 *     algorithm:      Cuba Cuhre (cubature)
 *     implementation: member function that expects x and y
 * 
 *  5) boundaries:     0<y0<3, -1<y1<1, -infinity<y2<3
 *     algorithm:      Cubature parallelised hcubature (cubature)
 *     implementation: lambda expression that expects x and y
 * 
 * The desired error limits as well as the maximal numbers of intervals and function evaluations, respectively, are chosen
 * such that the integral succeeds in cases 1, 3, 4 and 5, and fails in case 2. A corresponding warning will be written to
 * the standard output.
 */

double global_function_1 (const MultiDimInt::Arguments<2>& x, const MultiDimInt::Arguments<3>& y)	// function implementation of test case 1
{
	const double x0 = x[0];
	const double x1 = x[1];
	
	const double y0 = y[0];
	const double y1 = y[1];
	const double y2 = y[2];
	
	return x0*y0 * std::pow(x1-y1, 4.0) * std::exp(-y2*y2);
}

class FunctionObject2	// function implementation of test case 2
{
public:
	FunctionObject2 () {};
	
	double operator() (const MultiDimInt::Arguments<2>& x, const MultiDimInt::Arguments<3>& y) const
	{
		const double x0 = x[0];
		const double x1 = x[1];
		
		const double y0 = y[0];
		const double y1 = y[1];
		const double y2 = y[2];
		
		return x0*y0 * std::pow(x1-y1, 4.0) * std::exp(-y2*y2);
	};
};

class FunctionObject3	// function implementation of test case 3
{
public:
	FunctionObject3 (const MultiDimInt::Arguments<2>& x) : x0(x[0]), x1(x[1]) {};
	
	double operator() (const MultiDimInt::Arguments<3>& y) const
	{
		const double y0 = y[0];
		const double y1 = y[1];
		const double y2 = y[2];
		
		return x0*y0 * std::pow(x1-y1, 4.0) * std::exp(-y2*y2);
	};
	
private:
	double x0, x1;
};

class Class4	// function implementation of test case 4
{
public:
	Class4 () {};
	
	double member_function_4 (const MultiDimInt::Arguments<2>& x, const MultiDimInt::Arguments<3>& y) const
	{
		const double x0 = x[0];
		const double x1 = x[1];
		
		const double y0 = y[0];
		const double y1 = y[1];
		const double y2 = y[2];
		
		return x0*y0 * std::pow(x1-y1, 4.0) * std::exp(-y2*y2);
	};
};

int main()
{
	const MultiDimInt::Arguments<2> x = {2.0, 1.0};	// the fixed arguments x0=2 and x1=1
	
	const MultiDimInt::Arguments<3> lowerBounds2 = {0.0, -1.0, -2.0};							// the lower integral boundaries used in test case 2
	const MultiDimInt::Arguments<3> upperBounds2 = {3.0, 1.0, 3.0};								// the upper integral boundaries used in test case 2
	const MultiDimInt::Arguments<3> lowerBounds3 = {0.0, -1.0, -2.0};							// the lower integral boundaries used in test case 3
	const MultiDimInt::Arguments<3> upperBounds3 = {3.0, 1.0, MultiDimInt::PositiveInfinity};	// the upper integral boundaries used in test case 3
	const MultiDimInt::Arguments<3> lowerBounds4 = {0.0, -1.0, MultiDimInt::NegativeInfinity};	// the lower integral boundaries used in test case 4
	const MultiDimInt::Arguments<3> upperBounds4 = {3.0, 1.0, MultiDimInt::PositiveInfinity};	// the upper integral boundaries used in test case 4
	const MultiDimInt::Arguments<3> lowerBounds5 = {0.0, -1.0, MultiDimInt::NegativeInfinity};	// the lower integral boundaries used in test case 5
	const MultiDimInt::Arguments<3> upperBounds5 = {3.0, 1.0, 3.0};								// the upper integral boundaries used in test case 5
	
	const double absErr = 1e-10;			// absolute error limit
	const double relErr = 1e-4;				// relative error limit
	const std::size_t maxInterval = 1e2;	// maximal number of intervals allowed in the integration algorithm used in test case 1
	const int maxEval = 1e5;				// maximal number of function evaluations allowed in the integration algorithms used in test cases 2, 3, 4 and 5
	
	const MultiDimInt::GSLNestedCQUADAlgorithm alg1(absErr, relErr, maxInterval);			// integration algorithm used in test case 1
	const MultiDimInt::GSLMonteCarloVegasAlgorithm alg2(absErr, relErr, maxEval);			// integration algorithm used in test case 2
	const MultiDimInt::CubaDivonneAlgorithm alg3(absErr, relErr, maxEval);					// integration algorithm used in test case 3
	const MultiDimInt::CubaCuhreAlgorithm alg4(absErr, relErr, maxEval);					// integration algorithm used in test case 4
	const MultiDimInt::CubatureParallelHAdaptiveAlgorithm alg5(absErr, relErr, maxEval);	// integration algorithm used in test case 5
	
	const FunctionObject2 funcObj2;		// instance of 'FunctionObject2' used in test case 2
	const FunctionObject3 funcObj3(x);	// instance of 'FunctionObject3' used in test case 3 (gets initialized with fixed arguments x)
	const Class4 instance4;				// instance of 'Class4' used in test case 4
	const auto lambda5 = [](const MultiDimInt::Arguments<2> &x, const MultiDimInt::Arguments<3> &y)
	{
		const double x0 = x[0];
		const double x1 = x[1];

		const double y0 = y[0];
		const double y1 = y[1];
		const double y2 = y[2];

		return x0 * y0 * std::pow(x1 - y1, 4.0) * std::exp(-y2 * y2);
	}; // lambda expression used in test case 5

	const MultiDimInt::Integrator<2,3> integrator1(global_function_1, alg1, "Test case 1");						// integrator used in test case 1
	const MultiDimInt::Integrator<2,3> integrator2(funcObj2, alg2, "Test case 2");								// integrator used in test case 2
	const MultiDimInt::Integrator<0,3> integrator3(funcObj3, alg3, "Test case 3");								// integrator used in test case 3 (here, the first template parameter (number of fixed arguments) is 0, as x is not an argument of testFuncObj3 but was provided via its constructor)
	const MultiDimInt::Integrator<2,3> integrator4(&Class4::member_function_4, instance4, alg4, "Test case 4");	// integrator used in test case 4 (for member functions a different constructor syntax has to be used)
	const MultiDimInt::Integrator<2,3> integrator5(lambda5, alg5, "Test case 5");								// integrator used in test case 5
	
	double result1, result2, result3, result4, result5;	// the integration results will be written into these variables
	double error1, error2, error3, error4, error5;		// the respective errors will be written into these variables
	
	integrator1.integrate(x, result1, error1);								// performing the integration in test case 1 (if no integration boundaries are provided as arguments to the 'integrate' method, the integration is performed over the unit hypercube)
	integrator2.integrate(x, lowerBounds2, upperBounds2, result2, error2);	// performing the integration in test case 2
	integrator3.integrate(lowerBounds3, upperBounds3, result3, error3);		// performing the integration in test case 3 (here, 'integrate' expects no fixed arguments, as x is not an argument of testFuncObj3 but was provided via its constructor)
	integrator4.integrate(x, lowerBounds4, upperBounds4, result4, error4);	// performing the integration in test case 4
	integrator5.integrate(x, lowerBounds5, upperBounds5, result5, error5);	// performing the integration in test case 5
	
	std::cout << std::scientific	// write the integration results, their errors and the exact results to the standard output
			  << std::endl
			  << "Test case 1: The result of the numerical integration is " << result1 << " +- " << error1 << "  (exact result = " << std::sqrt(M_PI) * std::erf(1.0) / 10.0						<< ")." << std::endl
			  << "Test case 2: The result of the numerical integration is " << result2 << " +- " << error2 << "  (exact result = " << std::sqrt(M_PI) * (std::erf(2.0)+std::erf(3.0)) * 144.0/5.0	<< ")." << std::endl
			  << "Test case 3: The result of the numerical integration is " << result3 << " +- " << error3 << "  (exact result = " << std::sqrt(M_PI) * (std::erf(2.0)+1.0) * 144.0/5.0				<< ")." << std::endl
			  << "Test case 4: The result of the numerical integration is " << result4 << " +- " << error4 << "  (exact result = " << std::sqrt(M_PI) * 288.0/5.0									<< ")." << std::endl
			  << "Test case 5: The result of the numerical integration is " << result5 << " +- " << error5 << "  (exact result = " << std::sqrt(M_PI) * (1.0+std::erf(3.0)) * 144.0/5.0				<< ")." << std::endl
			  << std::endl;
	
	return 0;
}