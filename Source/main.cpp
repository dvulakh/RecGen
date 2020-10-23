
#include "poisson_pedigree.h"
#include "rec_gen.h"

#include <iostream>
using namespace std;

int main() {
	auto ped = poisson_pedigree(10, 3, 4, 4).build();
	//auto ped2 = poisson_pedigree::recover_dumped(ped->dump(), new poisson_pedigree());
	//cout << ped->dump() << endl << endl << ped->dump_extant();
	//cout << "\n=======================================\n\n" << ped2->dump();
	auto ped2 = poisson_pedigree::recover_dumped(ped->dump_extant(), new poisson_pedigree());
	auto rg = new rec_gen(ped2);
	rg->apply_rec_gen();
	return 0;
}
