#include "CubatureParallelAlgorithm.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

bool MultiDimInt::CubatureParallelAlgorithm::is_parallelized() const
{
	return true; // all vectorised Cubature integration algorithms are used to parallelize the sampling of the integrand
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

MultiDimInt::CubatureParallelAlgorithm::CubatureParallelAlgorithm(const double absErr, const double relErr, const std::size_t maxEval)
	: CubatureAlgorithm(absErr, relErr, maxEval)
{
}

int MultiDimInt::CubatureParallelAlgorithm::cubature_integrand(const unsigned dimInt, const std::size_t numPoints, const double *argsInt, void *cubatureData, const unsigned dimFunc, double *result)
{
	CubatureData *data = (CubatureData *)cubatureData;

#pragma omp parallel for schedule(static) // call the integrand for all points in parallel
	for (std::size_t i_point = 0; i_point < numPoints; ++i_point)
	{
		result[i_point] = data->Func(data->ArgsFix, &argsInt[i_point * dimInt]); // evaluate integrand
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private