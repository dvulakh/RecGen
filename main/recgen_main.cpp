
/********************************************************************
* Reads the extant population of a poisson pedigree from STDIN and
* writes a poisson pedigree rebuilt with REC-GEN to STDOUT.
********************************************************************/

#include "../source/poisson_pedigree.h"
#include "../source/rec_gen_recursive.h"
#include "../source/rec_gen_bp.h"
#include "../source/flags.h"

#include <iostream>

#define STOP_CHAR '~'

int main(int narg, char** args)
{

	// Construct pedigree from STDIN
	std::string extant_dump, line;
	std::getline(std::cin, extant_dump, STOP_CHAR);
	poisson_pedigree* ped = poisson_pedigree::recover_dumped(extant_dump, new poisson_pedigree());

	// Prepare the Rec-Gen object
	rec_gen* recrec = new rec_gen_recursive(ped);
	rec_gen* recgen = new rec_gen_quadratic(ped);
	rec_gen* recbas = new rec_gen_basic(ped);
	rec_gen* recbp  = new rec_gen_bp(ped);

	recrec->settings = recgen->settings = recbas->settings = LOG_WORK | LOG_DATA;

	// Flag definitions
	flag_reader fr;
	LOG_FLAG_READ(fr, recgen);
	fr.add_flag("sib", 'S', 1, [&](std::vector<std::string> v, void* p) {
		std::vector<double> parg;
		for (auto s : split_opts(v[0]))
			parg.push_back(std::stod(s));
		recgen->set_sib(parg);
	});
	fr.add_flag("cand", 'c', 1, [&](std::vector<std::string> v, void* p) {
		std::vector<double> parg;
		for (auto s : split_opts(v[0]))
			parg.push_back(std::stod(s));
		recgen->set_cand(parg);
	});
	fr.add_flag("rec", 'r', 1, [&](std::vector<std::string> v, void* p) { recgen->set_rec(std::stod(v[0])); });
	fr.add_flag("decay", 'y', 1, [&](std::vector<std::string> v, void* p) { recgen->set_dec(std::stod(v[0])); });
	fr.add_flag("richness", 'd', 1, [&](std::vector<std::string> v, void* p) { recgen->set_d(std::stoi(v[0])); });
	fr.add_flag("basic", 'B', 0, [&](std::vector<std::string> v, void* p) { recgen = recbas; });
	fr.add_flag("recursive", 'R', 0, [&](std::vector<std::string> v, void* p) { recgen = recrec; });
	fr.add_flag("bp", 'P', 0, [&](std::vector<std::string> v, void* p) { recgen = recbp; });
	fr.add_flag("epsilon", 'e', 1, [&](std::vector<std::string> v, void* p) { static_cast<rec_gen_bp*>(recgen)->set_epsilon(std::stod(v[0])); });
	fr.add_flag("memmode", 'm', 1, [&](std::vector<std::string> v, void* p) { static_cast<rec_gen_bp*>(recgen)->set_memory_mode(std::stoi(v[0])); });

	if (fr.read_flags(narg, args) != FLAGS_INPUT_SUCCESS) {
		std::cout << "Invalid commands" << std::endl;
		return 1;
	}

	// Run REC-GEN
	recgen->init()->apply_rec_gen();
	std::cout << recgen->get_pedigree()->dump() << std::endl;
	delete ped;
	return 0;

}
