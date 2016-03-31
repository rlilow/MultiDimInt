#include "CubaCuhreAlgorithm.h"

#include <string>

#include <cuba.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::CubaCuhreAlgorithm::CubaCuhreAlgorithm (const double absErr, const double relErr, const int maxEval) :
	CubaAlgorithm(absErr, relErr, maxEval),
	Key(0)
{}

MultiDimInt::CubaCuhreAlgorithm::CubaCuhreAlgorithm (const double absErr, const double relErr, const int maxEval,
													 const int flags, const int mineval, const std::string& statefile, void* spin,
													 const int key) :
	CubaAlgorithm(absErr, relErr, maxEval,
				  flags, mineval, statefile, spin),
	Key(key)
{}

MultiDimInt::CubaCuhreAlgorithm* MultiDimInt::CubaCuhreAlgorithm::clone () const
{
	return new CubaCuhreAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

bool MultiDimInt::CubaCuhreAlgorithm::cuba_integration (const int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const
{
	const int ncomp = 1;	// number of components of the integrand (always 1 in MultiDimInt::Function)
	const int nvec = 1;		// number of integration points sampled at the same time (always 1 in MultiDimInt::Function)
	
	int nregions;	// actual number of subregions needed (will not be used)
	int neval;		// actual number of integrand evaluations needed (will not be used)
	int fail;		// Cuba error code
	
	Cuhre(dimInt, ncomp,
		  &cuba_integrand, &cubaData, nvec,
		  RelErr, AbsErr,
		  Flags,
		  Mineval, MaxEval,
		  Key, Statefile.c_str(), Spin,
		  &nregions, &neval, &fail,
		  &value, &error, &prob);
	
	if ( fail == 0 )	// if integration succeeded, return 'true'
	{
		return true;
	}
	else				// if integration failed, return 'false'
	{
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private