
/********************************************************************
* Reads the extant population of a poisson pedigree from STDIN and
* writes information about its properties to STDOUT.
********************************************************************/

#include "../source/tree_analyze.h"

#include <iostream>
#include <sstream>

#define STOP_CHAR '~'

std::vector<std::string> split_opts(std::string s)
{
	std::stringstream sin(s);
	std::vector<std::string> opts;
	while (sin.good()) {
		std::string o;
		getline(sin, o, ',');
		opts.push_back(o);
	}
	return opts;
}

int main(int narg, char** args)
{

	// Construct pedigree from STDIN
	std::string extant_dump, line;
	std::getline(std::cin, extant_dump, STOP_CHAR);
	poisson_pedigree* ped = poisson_pedigree::recover_dumped(extant_dump, new poisson_pedigree());

	// Get analysis data through flags
	preprocess* prep = new preprocess(ped);
	flag_reader fr;
	fr.add_flag("badlca", 'L', 0, [&](std::vector<std::string> v, void* p) {
		auto bad_lca = bad_joint_LCAs(prep);
		for (int i = 1; i < ped->num_grade(); i++)
			std::cout << "Generation " << i << ":\t" << bad_lca[i].first << "/" << bad_lca[i].second << "\t" << (100 * bad_lca[i].first / std::max(1LL, bad_lca[i].second)) << "%\n";
		std::cout << std::endl;
	});
	fr.add_flag("blocks", 'B', 1, [&](std::vector<std::string> v, void* p) {
		auto opts = split_opts(v[0]);
		bool div = opts[0] != "0" && opts[0] != "";
		int gen = opts.size() > 1 && opts[1] != "" ? stoi(opts[1]) : 0;
		auto block_stat = block_share_stat(prep, gen);
		for (int i = 0; i < ped->num_grade(); i++)
			for (int j = 0; j < 3; j++) {
				for (int e : block_stat[i][j])
					std::cout << (div ? 100 * e / ped->num_blocks() : e) << " ";
				std::cout << std::endl;
			}
		std::cout << std::endl;
	});
	fr.add_flag("siblocks", 'b', 1, [&](std::vector<std::string> v, void* p) {
		bool div = v[0] != "0";
		auto block_stat = sib_block_share_stat(prep);
		for (auto& e : block_stat) {
			for (auto i : e)
				std::cout << (div ? 100 * i / ped->num_blocks() : i) << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	});
	if (fr.read_flags(narg, args) != FLAGS_INPUT_SUCCESS) {
		std::cout << "Invalid commands" << std::endl;
		return 1;
	}
	return 0;

}
