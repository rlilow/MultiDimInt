#include "CubaDivonneAlgorithm.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <cuba.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

MultiDimInt::CubaDivonneAlgorithm::CubaDivonneAlgorithm (const double absErr, const double relErr, const int maxEval) :
	CubaAlgorithm(absErr, relErr, maxEval),
	Seed(0),
	Key1(47),
	Key2(1),
	Key3(1),
	Maxpass(5),
	Border(0.0),
	Maxchisq(10.0),
	Mindeviation(0.25),
	Ngiven(0),
	Xgiven(),
	Nextra(0),
	Peakfinder(NULL)
{}

MultiDimInt::CubaDivonneAlgorithm::CubaDivonneAlgorithm (const double absErr, const double relErr, const int maxEval,
														 const int flags, const int mineval, const std::string& statefile, void* spin,
														 const int seed, const int key1, const int key2, const int key3, const int maxpass, const double border, const double maxchisq, const double mindeviation, const int ngiven, const std::vector<double>& xgiven, const int nextra, peakfinder_t peakfinder) :
	CubaAlgorithm(absErr, relErr, maxEval,
				  flags, mineval, statefile, spin),
	Seed(seed),
	Key1(key1),
	Key2(key2),
	Key3(key3),
	Maxpass(maxpass),
	Border(border),
	Maxchisq(maxchisq),
	Mindeviation(mindeviation),
	Ngiven(ngiven),
	Xgiven(xgiven),
	Nextra(nextra),
	Peakfinder(peakfinder)
{
	if ( Seed < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaDivonneAlgorithm Error: Value of 'seed' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Maxpass <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaDivonneAlgorithm Error: Value of 'maxpass' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Border < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaDivonneAlgorithm Error: Value of 'border' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Maxchisq <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaDivonneAlgorithm Error: Value of 'maxchisq' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Mindeviation <= 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaDivonneAlgorithm Error: Value of 'mindeviation' is not positive" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Ngiven < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaDivonneAlgorithm Error: Value of 'ngiven' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
	
	if ( Nextra < 0 )
	{
		std::cout << std::endl
				  << " MultiDimInt::CubaDivonneAlgorithm Error: Value of 'nextra' is negative" << std::endl
				  << std::endl;
		
		exit(EXIT_FAILURE);
	}
}

MultiDimInt::CubaDivonneAlgorithm* MultiDimInt::CubaDivonneAlgorithm::clone () const
{
	return new CubaDivonneAlgorithm(*this);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// protected

bool MultiDimInt::CubaDivonneAlgorithm::cuba_integration (const int dimInt, CubaData& cubaData, double& value, double& error, double& prob, std::string& furtherComment) const
{
	const int ncomp = 1;	// number of components of the integrand (always 1 in MultiDimInt::Function)
	const int nvec = 1;		// number of integration points sampled at the same time (always 1 in MultiDimInt::Function)
	
	const int ldxgiven = dimInt;	// offset between one point and the next in the array 'Xgiven' (always assumed to be given by 'dimInt', i.e. the dimension of a sample point)
	
	int nregions;	// actual number of subregions needed (will not be used)
	int neval;		// actual number of integrand evaluations needed (will not be used)
	int fail;		// Cuba error code
	
	Divonne(dimInt, ncomp,
			&cuba_integrand, &cubaData, nvec,
			RelErr, AbsErr,
			Flags, Seed,
			Mineval, MaxEval,
			Key1, Key2, Key3,
			Maxpass, Border,
			Maxchisq, Mindeviation,
			Ngiven, ldxgiven, const_cast<double*>(Xgiven.data()),
			Nextra, Peakfinder,
			Statefile.c_str(), Spin,
			&nregions, &neval, &fail,
			&value, &error, &prob);
	
	if ( fail == 0 )	// if integration succeeded, return 'true'
	{
		return true;
	}
	else				// if integration failed, write the estimated increase of integrand evaluations 'maxEval' needed to reach the specified tolerance (contained in 'fail') into 'furtherComment' and return 'false'
	{
		furtherComment = std::string("\n")
					   + std::string("	-Estimated increase of integrand evaluations needed to succeed: ") + std::to_string(fail);
		
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// private