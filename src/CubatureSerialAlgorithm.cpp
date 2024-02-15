#include "CubatureSerialAlgorithm.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

bool MultiDimInt::CubatureSerialAlgorithm::is_parallelized() const
{
	return false; // all non-vectorised Cubature integration algorithms use only a single thread
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

MultiDimInt::CubatureSerialAlgorithm::CubatureSerialAlgorithm(const double absErr, const double relErr, const std::size_t maxEval)
	: CubatureAlgorithm(absErr, relErr, maxEval)
{
}

int MultiDimInt::CubatureSerialAlgorithm::cubature_integrand(const unsigned dimInt, const double *argsInt, void *cubatureData, const unsigned dimFunc, double *result)
{
	CubatureData *data = (CubatureData *)cubatureData;

	result[0] = data->Func(data->ArgsFix, argsInt); // evaluate integrand

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private