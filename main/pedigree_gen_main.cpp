
/********************************************************************
* Stochastically generates a poisson pedigree based on properties
* read from command-line arguments and prints it to STDOUT.
********************************************************************/

#include "../source/poisson_pedigree.h"
#include "../source/flags.h"

#include <iostream>

#define STOP_CHAR '~'

int main(int narg, char** args)
{

	// Read pedigree properties
	std::string arg;
	for (int i = 1; i < narg; i++)
		arg += std::string(args[i]) + " ";
	poisson_pedigree* ped = poisson_pedigree::recover_dumped(arg, new poisson_pedigree());

	// Generate and print pedigree
	std::cout << ped->build()->dump_extant() << std::endl << STOP_CHAR << std::endl << ped->dump() << std::endl;
	return 0;

}
