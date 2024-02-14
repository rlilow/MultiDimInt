#include "CubaSuaveAlgorithm.hpp"

#include <iostream>
#include <string>

#include <cuba.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::CubaSuaveAlgorithm::CubaSuaveAlgorithm (const double absErr, const double relErr, const int maxEval) :
	CubaAlgorithm(absErr, relErr, maxEval),
	Seed(0),
	Nnew(1000),
	Nmin(2),
	Flatness(50.0)
{}

MultiDimInt::CubaSuaveAlgorithm::CubaSuaveAlgorithm (const double absErr, const double relErr, const int maxEval,
													 const int flags, const int mineval, const std::string& statefile, void* spin,
													 const int seed, const int nnew, const int nmin, const double flatness) :
	CubaAlgorithm(absErr, relErr, maxEval,
				  flags, mineval, statefile, spin),
	Seed(seed),
	Nnew(nnew),
	Nmin(nmin),
	Flatness(flatness)
{
	if ( Seed < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaSuaveAlgorithm Error: Value of 'seed' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Nnew <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaSuaveAlgorithm Error: Value of 'nnew' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Nmin <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaSuaveAlgorithm Error: Value of 'nmin' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Flatness <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaSuaveAlgorithm Error: Value of 'flatness' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::CubaSuaveAlgorithm* MultiDimInt::CubaSuaveAlgorithm::clone () const
{
	return new CubaSuaveAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

bool MultiDimInt::CubaSuaveAlgorithm::cuba_integration (const int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const
{
	const int ncomp = 1;	// number of components of the integrand (always 1 in MultiDimInt::Function)
	const int nvec = 1;		// number of integration points sampled at the same time (always 1 in MultiDimInt::Function)
	
	int nregions;	// actual number of subregions needed (will not be used)
	int neval;		// actual number of integrand evaluations needed (will not be used)
	int fail;		// Cuba error code
	
	Suave(dimInt, ncomp,
		  &cuba_integrand, &cubaData, nvec,
		  RelErr, AbsErr,
		  Flags, Seed,
		  Mineval, MaxEval,
		  Nnew, Nmin,
		  Flatness, Statefile.c_str(), Spin,
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