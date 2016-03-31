#include "../src/BindMemberFunction.h"

#include <cmath>
#include <functional>
#include <iostream>

/**
 * BindMemberFunction demo:
 * 
 * Binding a member function to an instance of its host class.
 * 
 * Note that this is not related to the actual use of MultiDimInt. Instead, it demonstrates the use of an auxiliary function
 * which is used internally in this library. This function, however, may be useful on its own.
 * 
 * If one wants to use a member function in a situation where a normal function is expected, the wrapper bind_member_function_to_object
 * can be used. This wrapper only expects a pointer to the member function and a reference to an instance of its host class.
 * 
 * Here, this is demonstrated for a simple test class whose member function is given by f(x,y) = x + y^n for some power
 * n specified in the constructor. For the two cases n=1 and n=2 the member function f(6,3) gets evaluated once through
 * an instance of its host class and once by calling the function obtained from binding it to that instance. As expected,
 * both approaches yield the same result.
 */

class Class	// simple test class
{
public:
	Class (int n) : N(n) {};
	
	double member_function (double x, double y) const
	{
		return x + std::pow(y, N);
	};
	
private:
	int N;
};

int main()
{
	const Class instance1(1);	// instance of 'Class' with n=1
	const Class instance2(2);	// instance of 'Class' with n=2
	
	const auto function_1 = MultiDimInt::bind_member_function_to_object(&Class::member_function, instance1);	// binding the member function to an instance of its host class with n=1 ('auto' automatically deduces the type of function_1)
	const auto function_2 = MultiDimInt::bind_member_function_to_object(&Class::member_function, instance2);	// binding the member function to an instance of its host class with n=1 ('auto' automatically deduces the type of function_2)
	
	const double x = 6.0;	// x argument of interest
	const double y = 3.0;	// y argument of interest
	
	std::cout << std::endl
			  << "Usual way to call member function: instance1.member_function(x,y) = "	<< instance1.member_function(x, y)	<< std::endl
			  << "Calling the bound member function: function_1(x,y) = "				<< function_1(x,y)					<< std::endl
			  << std::endl
			  << "Usual way to call member function: instance2.member_function(x,y) = "	<< instance2.member_function(x, y)	<< std::endl
			  << "Calling the bound member function: function_2(x,y) = "				<< function_2(x,y)					<< std::endl
			  << std::endl;
	
	return 0;
}