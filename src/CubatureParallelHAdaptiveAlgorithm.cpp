#include "CubatureParallelHAdaptiveAlgorithm.hpp"

#include <vector>

#include <cubature.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::CubatureParallelHAdaptiveAlgorithm::CubatureParallelHAdaptiveAlgorithm(const double absErr, const double relErr, const std::size_t maxEval)
	: CubatureParallelAlgorithm(absErr, relErr, maxEval)
{
}

MultiDimInt::CubatureParallelHAdaptiveAlgorithm *MultiDimInt::CubatureParallelHAdaptiveAlgorithm::clone() const
{
	return new CubatureParallelHAdaptiveAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

bool MultiDimInt::CubatureParallelHAdaptiveAlgorithm::cubature_integration(const unsigned dimInt, CubatureData &cubatureData, double &value, double &error, std::string &furtherComment) const
{
	const unsigned dimFunc = 1; // number of components of the integrand (always 1 in MultiDimInt::Function)

	std::vector<double> lowerBound(dimInt, 0.0); // this must not be const, as hcubature_v expects the integration boundaries as arrays of non-const double values
	std::vector<double> upperBound(dimInt, 1.0);

	int fail = hcubature_v(dimFunc, &cubature_integrand, &cubatureData,
						   dimInt, lowerBound.data(), upperBound.data(),
						   MaxEval, AbsErr, RelErr,
						   ERROR_INDIVIDUAL, // placeholder value, since the error norm argument is ignored for integrands with only 1 component
						   &value, &error);

	if (fail == 0) // if integration succeeded, return 'true'
	{
		return true;
	}
	else // if integration failed, return 'false'
	{
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private