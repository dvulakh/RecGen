
#include "poisson_pedigree.h"
#include "rec_gen_basic.h"
#include "tree_diff_basic.h"

#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;

int main() {
	// auto ped = poisson_pedigree(50, 3, 4, 4).build();
	// //auto ped2 = poisson_pedigree::recover_dumped(ped->dump(), new poisson_pedigree());
	// //cout << ped->dump() << endl << endl << ped->dump_extant();
	// //cout << "\n=======================================\n\n" << ped2->dump();
	//*
	ifstream fin ("good_extant");
	ifstream fin2("good_dump");
	string s;
	getline (fin, s, (char)fin.eof());
	auto ped2 = poisson_pedigree::recover_dumped(s, new poisson_pedigree());
	getline (fin2, s, (char)fin2.eof());
	auto ped = poisson_pedigree::recover_dumped(s, new poisson_pedigree());
	//*/
	/*
	//auto ped = poisson_pedigree(100, 6, 3, 15).build();
	auto ped = poisson_pedigree(10, 5, 2, 4).build();
	cout << ped->dump() << endl;
	cout << ped->dump_extant() << endl;
	auto ped2 = poisson_pedigree::recover_dumped(ped->dump_extant(), new poisson_pedigree());
	//*/
	//*
	auto rg = new rec_gen_basic(ped2);
	cout << ped2->dump_extant() << endl;
	rg->apply_rec_gen();
	cout << endl << rg->get_pedigree()->dump() << endl;
	auto diff = (new tree_diff_basic(ped, ped2))->topology_biject()->blocks_check();
	printf("\nNodes correct:    %d/%d \t(%d%%)\nEdges correct:    %d/%d\t(%d%%)\nBlocks attempted: %d/%d\t(%d%%)\nBlocks correct:   %d/%d\t(%d%%/%d%%)\n",
		diff->nodes_correct, diff->nodes_total, 100 * diff->nodes_correct / diff->nodes_total,
		diff->edges_correct, diff->edges_total, 100 * diff->edges_correct / diff->edges_total,
		diff->blocks_attempted, diff->blocks_total, 100 * diff->blocks_attempted / diff->blocks_total,
		diff->blocks_correct, diff->blocks_total, 100 * diff->blocks_correct / diff->blocks_total, 100 * diff->blocks_correct / diff->blocks_attempted);
	//*/
	return 0;
}
