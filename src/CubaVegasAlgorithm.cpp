#include "CubaVegasAlgorithm.hpp"

#include <cmath>
#include <iostream>
#include <string>

#include <cuba.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::CubaVegasAlgorithm::CubaVegasAlgorithm (const double absErr, const double relErr, const int maxEval) :
	CubaAlgorithm(absErr, relErr, maxEval),
	Seed(0),
	Nstart(1000),
	Nincrease(500),
	Nbatch(1000),
	Gridno(0)
{}

MultiDimInt::CubaVegasAlgorithm::CubaVegasAlgorithm (const double absErr, const double relErr, const int maxEval,
													 const int flags, const int mineval, const std::string& statefile, void* spin,
													 const int seed, const int nstart, const int nincrease, const int nbatch, const int gridno) :
	CubaAlgorithm(absErr, relErr, maxEval,
				  flags, mineval, statefile, spin),
	Seed(seed),
	Nstart(nstart),
	Nincrease(nincrease),
	Nbatch(nbatch),
	Gridno(gridno)
{
	if ( Seed < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaSuaveAlgorithm Error: Value of 'seed' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Nstart <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaVegasAlgorithm Error: Value of 'nstart' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Nincrease <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaVegasAlgorithm Error: Value of 'nincrease' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Nbatch <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaVegasAlgorithm Error: Value of 'nbatch' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( std::abs(Gridno) > 10 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaVegasAlgorithm Error: Invalid value of 'gridno'" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::CubaVegasAlgorithm* MultiDimInt::CubaVegasAlgorithm::clone () const
{
	return new CubaVegasAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

bool MultiDimInt::CubaVegasAlgorithm::cuba_integration (const int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const
{
	const int ncomp = 1;	// number of components of the integrand (always 1 in MultiDimInt::Function)
	const int nvec = 1;		// number of integration points sampled at the same time (always 1 in MultiDimInt::Function)
	
	int neval;	// actual number of integrand evaluations needed (will not be used)
	int fail;	// Cuba error code
	
	Vegas(dimInt, ncomp,
		  &cuba_integrand, &cubaData, nvec,
		  RelErr, AbsErr,
		  Flags, Seed,
		  Mineval, MaxEval,
		  Nstart, Nincrease, Nbatch,
		  Gridno, Statefile.c_str(), Spin,
		  &neval, &fail,
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