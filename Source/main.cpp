
#include "poisson_pedigree.h"
#include "rec_gen_basic.h"

#include <iostream>
#include <fstream>
using namespace std;

int main() {
	// auto ped = poisson_pedigree(50, 3, 4, 4).build();
	// //auto ped2 = poisson_pedigree::recover_dumped(ped->dump(), new poisson_pedigree());
	// //cout << ped->dump() << endl << endl << ped->dump_extant();
	// //cout << "\n=======================================\n\n" << ped2->dump();
	/*
	ifstream fin ("good_extant");
	string s;
	getline (fin, s, (char)fin.eof());
	auto ped2 = poisson_pedigree::recover_dumped(s, new poisson_pedigree());
	//*/
	//*
	auto ped = poisson_pedigree(100, 5, 2, 9).build();
	auto ped2 = poisson_pedigree::recover_dumped(ped->dump_extant(), new poisson_pedigree());
	cout << ped->dump() << endl;
	//*/
	auto rg = new rec_gen_basic(ped2);
	cout << ped2->dump_extant() << endl;
	rg->apply_rec_gen();
	cout << rg->get_pedigree()->dump() << endl;
	return 0;
}
