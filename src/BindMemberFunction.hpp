#ifndef MULTIDIMINT_BIND_MEMBER_FUNCTION_H
#define MULTIDIMINT_BIND_MEMBER_FUNCTION_H

#include <functional>

namespace MultiDimInt
{
	/**
	  * Wrapper that turns a general member function into a \c std::function such that it can be used in any situation
	  * where a usual function is expected. It takes a pointer to a member function \a memberFunctionPointer as well
	  * as a reference to an instance \a object of its host class \a Class as arguments.
	  * 
	  * There are no restrictions to the number and types of the arguments, \a ArgumentTypes, or the type of the return
	  * value, \a ReturnType, of the member function.
	  * 
	  * Author: Robert Lilow, ITA, ZAH, Heidelberg University (2016)
	  */
	template <class Class, typename ReturnType, typename... ArgumentTypes>
	std::function<ReturnType(ArgumentTypes...)> bind_member_function_to_object(ReturnType(Class::* const memberFunctionPointer)(ArgumentTypes...), Class& object)
	{
		Class* const objectPointer = &object;
		
		return [=] (ArgumentTypes... arguments) {return (objectPointer->*memberFunctionPointer)(arguments...);};
	}
	
	/**
	  * Wrapper that turns a general \c const member function into a \c std::function such that it can be used in any
	  * situation where a usual function is expected. It takes a pointer to a \c const member function \a constMemberFunctionPointer
	  * as well as a reference to a \c const instance \a constObject of its host class \a Class as arguments.
	  * 
	  * There are no restrictions to the number and types of the arguments, \a ArgumentTypes, or the type of the return
	  * value, \a ReturnType, of the \c const member function.
	  * 
	  * Author: Robert Lilow, ITA, ZAH, Heidelberg University (2016)
	  */
	template <class Class, typename ReturnType, typename... ArgumentTypes>
	std::function<ReturnType(ArgumentTypes...)> bind_member_function_to_object(ReturnType(Class::* const constMemberFunctionPointer)(ArgumentTypes...) const, const Class& constObject)
	{
		const Class* const constObjectPointer = &constObject;
		
		return [=] (ArgumentTypes... arguments) {return (constObjectPointer->*constMemberFunctionPointer)(arguments...);};
	}
}

#endif