#include "CubatureAlgorithm.hpp"

#include <climits>
#include <iostream>
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::Algorithm::Result MultiDimInt::CubatureAlgorithm::run(const InternalIntegrand &func, const std::size_t dimInt, const double *argsFix) const
{
	if (dimInt > INT_MAX) // Cubature algorithms only accept an 'unsigned' as the number of integration variables, not a potentially larger 'size_t'; in practice 'dimInt' will of course never exceed the largest possible 'unsigned' value, but explicitly checking this won't hurt
	{
		std::cout << std::endl
				  << " MultiDimInt::CubatureAlgorithm Error: Number of integration variables to large" << std::endl
				  << std::endl;

		exit(EXIT_FAILURE);
	}

	Algorithm::Result integral = {false, 0.0, 0.0, ""};

	CubatureData cubatureData(func, argsFix);

	bool integrationSucceeded;		// if the integration succeeds, this is set to 'true', otherwise to 'false'
	std::string furtherComment(""); // if the integration fails, an additional comment may be written to this string

	integrationSucceeded = cubature_integration(static_cast<int>(dimInt), cubatureData, integral.Value, integral.Error, furtherComment);

	if (not integrationSucceeded) // if integration failed, write Cubature error message and 'furtherComment' into 'integral.Comment'
	{
		integral.Failed = true;

		integral.Comment = std::string("	-Cubature error: Failed to reach the specified tolerance") + std::string("\n") + furtherComment;
	}

	return integral;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

MultiDimInt::CubatureAlgorithm::CubatureAlgorithm(const double absErr, const double relErr, const std::size_t maxEval)
	: Algorithm(absErr, relErr),
	  MaxEval(maxEval)
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private