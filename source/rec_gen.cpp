
/********************************************************************
* Implements the general structure of the Rec-Gen algorithm.
* Other implementation details are filled in by child classes.
********************************************************************/

#include "rec_gen.h"

// Parameter defaults
#define DEFAULT_SIB 0.21
#define DEFAULT_REC 0.99
#define DEFAULT_DEC 0.85
#define DEFAULT_D 3

// The Rec-Gen algorithm (Algorithm 1)
poisson_pedigree* rec_gen::apply_rec_gen()
{
	/// Reset the pedigree to the extant population
	start_time = std::chrono::high_resolution_clock::now();
	WPRINT(PRINT_HEADER("REC-GEN BEGINS"))
	ped->reset();
	/// Rebuild each grade
	while (!ped->done()) {
		/// Build the next generation
		WPRINT(PRINT_HEADER("NEW GENERATION"))
		WPRINT("Conducting siblinghood test")
		hypergraph* G = test_siblinghood();
		WPRINT("Assigning parents")
		assign_parents(G);
		delete G;
		/// Gather genetic information
		for (coupled_node* v : *ped) {
			WPRINTF("Collecting symbols for couple %lld", v->get_id())
			collect_symbols(v);
		}
	}
	WPRINT(PRINT_HEADER("DONE"))
	/// Set the founders as their own parents
	for (coupled_node *v : *ped)
		(*v)[0]->assign_par(v), (*v)[1]->assign_par(v);
	/// Return the pedigree
	return ped;
}

// Initialization and construction of rec-gen object
/// Initialize given all info
void rec_gen::init(poisson_pedigree* ped, std::string work_log, std::string data_log, double sib, double cand, double decay, double rec, int d, long long settings)
{
	this->ped = ped;
	this->work_path = work_log;
	this->data_path = data_log;
	this->work_log = NULL;
	this->data_log = NULL;
	this->sib = sib;
	this->cand = cand;
	this->decay = decay;
	this->rec = rec;
	this->d = d;
	this->settings = settings;
	this->init();
}
/// Initialize based on current members
rec_gen* rec_gen::init()
{
	/// Open files
	if (IS(LOG_WORK))
		delete this->work_log, this->work_log = std::fopen(this->work_path.c_str(), "w");
	if (IS(LOG_DATA))
		delete this->data_log, this->data_log = std::fopen(this->data_path.c_str(), "w");
	return this;
}
/// Construct given all info
rec_gen::rec_gen(poisson_pedigree* ped, std::string work_log, std::string data_log, double sib, double cand, double decay, double rec, int d, long long settings)
{ init(ped, work_log, data_log, sib, cand, decay, rec, d, settings); }
/// Construct given pedigree
rec_gen::rec_gen(poisson_pedigree* ped)
{ init(ped, "logs/rec-gen.log", "logs/rec-gen.dat", DEFAULT_SIB, DEFAULT_SIB, DEFAULT_REC, DEFAULT_DEC, DEFAULT_D, 0); }

// Access pedigree
poisson_pedigree* rec_gen::get_pedigree() { return this->ped; }

// Mutators
rec_gen* rec_gen::set_cand(double cand) { this->cand = cand; return this; }
rec_gen* rec_gen::set_sib(double sib) { this->sib = sib; return this; }
rec_gen* rec_gen::set_rec(double rec) { this->rec = rec; return this; }
rec_gen* rec_gen::set_dec(double decay) { this->decay = decay; return this; }
rec_gen* rec_gen::set_d(int d) { this->d = d; return this; }
