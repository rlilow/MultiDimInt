#include "../MultiDimInt.h"

#include <cmath>
#include <iostream>

/**
 * Nesting_different_algorithms demo:
 * 
 * Integral of the function f(x;y) = f(x0,x1;y0,y1,y2) = x0*y0 * (x1-y1)^4 * exp(-y2^2) over its arguments y0, y1 and y2,
 * each from 0 to 1, for fixed argument x0=2 and x1=1, using nesting of different integration algorithms.
 * 
 * To get an idea of the different possibilities, 3 exemplary cases are demonstrated:
 * 
 *  1) algorithm:      Cuba Cuhre for y0, y1 and y2
 *     implementation: global function that expects x as a fixed and y as integration arguments
 * 
 *  2) algorithm:      Cuba Cuhre for y0 and y1, and GSL CQUAD for y2
 *     implementation: 2 nested global functions where the 1st expects x as fixed and y0 as well as y1 as integration arguments, and the 2nd expects x, y0 as well as y1 as fixed and y2 as integration arguments
 * 
 *  3) algorithm:      GSL CQUAD for y0, GSL QAG for y1 and GSL QNG for y2
 *     implementation: 3 nested global functions where the 1st expects x as fixed and y0 as integration arguments, the 2nd expects x as well as y0 as fixed and y1 as integration arguments, and the 3rd expects x, y0 as well as y1 as fixed and y2 as integration arguments
 * 
 * Future versions of MultiDimInt may allow for a more direct way to handle the nesting of different algorithms.
 */

const double absErr = 1e-10;			// absolute error limit
const double relErr = 1e-4;				// relative error limit
const std::size_t maxInterval = 1e2;	// maximal number of intervals allowed in the GSL CQUAD and GSL QAG integration algorithms
const int maxEval = 1e5;				// maximal number of function evaluations allowed in the Cuba Cuhre integration algorithm

const MultiDimInt::CubaCuhreAlgorithm alg1_1(absErr, relErr, maxEval);			// integration algorithm used in test case 1
const MultiDimInt::CubaCuhreAlgorithm alg2_1(absErr, relErr, maxEval);			// integration algorithm used to integrate the 1st function in test case 2
const MultiDimInt::GSLNestedCQUADAlgorithm alg2_2(absErr, relErr, maxInterval);	// integration algorithm used to integrate the 2nd function in test case 2
const MultiDimInt::GSLNestedCQUADAlgorithm alg3_1(absErr, relErr, maxInterval);	// integration algorithm used to integrate the 1st function in test case 3
const MultiDimInt::GSLNestedQAGAlgorithm alg3_2(absErr, relErr, maxInterval);	// integration algorithm used to integrate the 2nd function in test case 3
const MultiDimInt::GSLNestedQNGAlgorithm alg3_3(absErr, relErr);				// integration algorithm used to integrate the 3rd function in test case 3

double function1_1 (const MultiDimInt::Arguments<2>& fixedArgs_x, const MultiDimInt::Arguments<3>& intArgs_y0_y1_y2)	// function for test case 1
{
	const double x0 = fixedArgs_x[0];
	const double x1 = fixedArgs_x[1];
	
	const double y0 = intArgs_y0_y1_y2[0];
	const double y1 = intArgs_y0_y1_y2[1];
	const double y2 = intArgs_y0_y1_y2[2];
	
	return x0*y0 * std::pow(x1-y1, 4.0) * std::exp(-y2*y2);	// evaluate f(x,y)
}

double function2_2 (const MultiDimInt::Arguments<4>& fixedArgs_x_y0_y1, const MultiDimInt::Arguments<1>& intArgs_y2)	// 2nd function for test case 2
{
	const double x0 = fixedArgs_x_y0_y1[0];
	const double x1 = fixedArgs_x_y0_y1[1];
	const double y0 = fixedArgs_x_y0_y1[2];
	const double y1 = fixedArgs_x_y0_y1[3];
	
	const double y2 = intArgs_y2[0];
	
	return x0*y0 * std::pow(x1-y1, 4.0) * std::exp(-y2*y2);	// evaluate f(x,y)
}

double function2_1 (const MultiDimInt::Arguments<2>& fixedArgs_x, const MultiDimInt::Arguments<2>& intArgs_y0_y1)	// 1st function for test case 2
{
	const double x0 = fixedArgs_x[0];
	const double x1 = fixedArgs_x[1];
	
	const double y0 = intArgs_y0_y1[0];
	const double y1 = intArgs_y0_y1[1];
	
	const MultiDimInt::Arguments<4> fixedArgs_x_y0_y1 = {x0, x1, y0, y1};	// fixed arguments expected by the 2nd function in test case 2
	
	const MultiDimInt::Integrator<4,1> integrator2_2(function2_2, alg2_2, "integrating function2_2");	// integrator used to integrate the 2nd function in test case 2
	
	double result2_2;	// the result of integrating the 2nd function in test case 2 will be written into this variable
	double error2_2;	// the respective error will be written into this variable (will be discarded)
	
	integrator2_2.integrate(fixedArgs_x_y0_y1, result2_2, error2_2);	// performing the integration of the 2nd function in test case 2
	
	return result2_2;
}

double function3_3 (const MultiDimInt::Arguments<4>& fixedArgs_x_y0_y1, const MultiDimInt::Arguments<1>& intArgs_y2)	// 3rd function for test case 3
{
	const double x0 = fixedArgs_x_y0_y1[0];
	const double x1 = fixedArgs_x_y0_y1[1];
	const double y0 = fixedArgs_x_y0_y1[2];
	const double y1 = fixedArgs_x_y0_y1[3];
	
	const double y2 = intArgs_y2[0];
	
	return x0*y0 * std::pow(x1-y1, 4.0) * std::exp(-y2*y2);	// evaluate f(x,y)
}

double function3_2 (const MultiDimInt::Arguments<3>& fixedArgs_x_y0, const MultiDimInt::Arguments<1>& intArgs_y1)	// 2nd function for test case 3
{
	const double x0 = fixedArgs_x_y0[0];
	const double x1 = fixedArgs_x_y0[1];
	const double y0 = fixedArgs_x_y0[2];
	
	const double y1 = intArgs_y1[0];
	
	const MultiDimInt::Arguments<4> fixedArgs_x_y0_y1 = {x0, x1, y0, y1};	// fixed arguments expected by the 3rd function in test case 3
	
	const MultiDimInt::Integrator<4,1> integrator3_3(function3_3, alg3_3, "integrating function3_3");	// integrator used to integrate the 3rd function in test case 3
	
	double result3_3;	// the result of integrating the 3rd function in test case 3 will be written into this variable
	double error3_3;	// the respective error will be written into this variable (will be discarded)
	
	integrator3_3.integrate(fixedArgs_x_y0_y1, result3_3, error3_3);	// performing the integration of the 3rd function in test case 3
	
	return result3_3;
}

double function3_1 (const MultiDimInt::Arguments<2>& fixedArgs_x, const MultiDimInt::Arguments<1>& intArgs_y0)	// 1st function for test case 3
{
	const double x0 = fixedArgs_x[0];
	const double x1 = fixedArgs_x[1];
	
	const double y0 = intArgs_y0[0];
	
	const MultiDimInt::Arguments<3> fixedArgs_x_y0 = {x0, x1, y0};	// fixed arguments expected by the 2nd function in test case 3
	
	const MultiDimInt::Integrator<3,1> integrator3_2(function3_2, alg3_2, "integrating function3_2");	// integrator used to integrate the 2nd function in test case 3
	
	double result3_2;	// the result of integrating the 2nd function in test case 3 will be written into this variable
	double error3_2;	// the respective error will be written into this variable (will be discarded)
	
	integrator3_2.integrate(fixedArgs_x_y0, result3_2, error3_2);	// performing the integration of the 2nd function in test case 3
	
	return result3_2;
}

int main()
{
	const MultiDimInt::Arguments<2> x = {2.0, 1.0};	// the fixed arguments x0=2 and x1=1
	
	const MultiDimInt::Integrator<2,3> integrator1_1(function1_1, alg1_1, "integrating function1_1");	// integrator used to integrate the function in test case 1
	const MultiDimInt::Integrator<2,2> integrator2_1(function2_1, alg2_1, "integrating function2_1");	// integrator used to integrate the 1st function in test case 2
	const MultiDimInt::Integrator<2,1> integrator3_1(function3_1, alg3_1, "integrating function3_1");	// integrator used to integrate the 1st function in test case 3
	
	double result1_1, result2_1, result3_1;	// the integration results will be written into these variables
	double error1_1, error2_1, error3_1;	// the respective errors will be written into these variables
	
	integrator1_1.integrate(x, result1_1, error1_1);	// performing the integration of the function in test case 1
	integrator2_1.integrate(x, result2_1, error2_1);	// performing the integration of the 1st function in test case 2
	integrator3_1.integrate(x, result3_1, error3_1);	// performing the integration of the 1st function in test case 3
	
	std::cout << std::scientific	// write the integration results, their errors and the exact result to the standard output
			  << std::endl
			  << "Test case 1: The result of the numerical integration is " << result1_1 << " +- " << error1_1			<< "." << std::endl
			  << "Test case 2: The result of the numerical integration is " << result2_1 << " +- " << error2_1			<< "." << std::endl
			  << "Test case 3: The result of the numerical integration is " << result3_1 << " +- " << error3_1			<< "." << std::endl
			  << "For comparison, the exact result is                     " << std::sqrt(M_PI) * std::erf(1.0) / 10.0	<< "." << std::endl
			  << std::endl;
	
	return 0;
}