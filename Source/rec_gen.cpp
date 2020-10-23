
/********************************************************************
* Implements the general structure of the Rec-Gen algorithm.
* Other implementation details are filled in by child classes.
********************************************************************/

#include "rec_gen.h"

// The Rec-Gen algorithm (Algorithm 1)
poisson_pedigree* rec_gen::apply_rec_gen()
{
	/// Reset the pedigree to the extant population
	start_time = std::chrono::high_resolution_clock::now();
	WPRINT(PRINT_HEADER("BEGIN SEARCHING"))
	ped->reset();
	/// Rebuild each grade
	while (!ped->done()) {
		/// Build the next generation
		WPRINT(PRINT_HEADER("NEW GENERATION"))
		WPRINT("Conducting siblinghood test")
		hypergraph* G = test_siblinghood();
		WPRINT("Assigning parents")
		assign_parents(G);
		/// Gather genetic information
		for (coupled_node* v : *ped) {
			WPRINTF("Collecting symbols for couple %lld", v->get_id())
			collect_symbols(v);
		}
	}
	WPRINT(PRINT_HEADER("DONE"))
	/// Return the pedigree
	return ped;
}

// Initialization and construction of rec-gen object
/// Initialize given all info
void rec_gen::init(poisson_pedigree* ped, std::string work_log, std::string data_log, long long settings)
{
	this->ped = ped;
	if (IS(LOG_WORK))
		this->work_log = std::fopen(work_log.c_str(), "w");
	if (IS(LOG_DATA))
		this->data_log = std::fopen(data_log.c_str(), "W");
	this->settings = settings;
}
/// Construct given all info
rec_gen::rec_gen(poisson_pedigree* ped, std::string work_log, std::string data_log, long long settings)
{ init(ped, work_log, data_log, settings); }
/// Construct given pedigree
rec_gen::rec_gen(poisson_pedigree* ped)
{ init(ped, "rec-gen.log", "rec-gen.dat", ~0); }

// Access pedigree
poisson_pedigree* rec_gen::get_pedigree() { return this->ped; }
