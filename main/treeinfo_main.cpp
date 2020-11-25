
/********************************************************************
* Reads the extant population of a poisson pedigree from STDIN and
* writes information about its properties to STDOUT.
********************************************************************/

#include "../source/tree_analyze.h"

#include <iostream>

#define STOP_CHAR '~'

int main(int narg, char** args)
{

	// Construct pedigree from STDIN
	std::string extant_dump, line;
	std::getline(std::cin, extant_dump, STOP_CHAR);
	poisson_pedigree* ped = poisson_pedigree::recover_dumped(extant_dump, new poisson_pedigree());

	// Get analysis data
	preprocess* prep = new preprocess(ped);
	auto bad_lca = bad_joint_LCAs(prep);

	// Output
	for (int i = 1; i < ped->num_grade(); i++)
		std::cout << "Generation " << i << ":\t" << bad_lca[i].first << "/" << bad_lca[i].second << "\t" << (100 * bad_lca[i].first / std::max(1LL, bad_lca[i].second)) << "%\n";
	return 0;

}
