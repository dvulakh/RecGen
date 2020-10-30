
/********************************************************************
* Defines the general structure of the Rec-Gen algorithm.
* Implementation details are filled in by child classes.
********************************************************************/

#ifndef REC_GEN_H
#define REC_GEN_H

#include "poisson_pedigree.h"

#include <chrono>
#include <cstdio>
#include <set>

// Accessing the settings bits of a rec_gen object
#define IS(S) (settings & (S))
#define LOG_WORK 1LL << 0
#define LOG_DATA 1LL << 1
#define VER_WORK 1LL << 2
#define VER_DATA 1LL << 3

// Logging macros
#define TPLUS(t0) (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - (t0)) / 1000000000.0)
#define PRINT_HEADER(s) (std::string(30, '=') + " " + s + " " + std::string(30, '='))
#define GPRINTF(s, ...) (std::string("[%f]\t") + s + std::string("\n")).c_str(), TPLUS(this->start_time), __VA_ARGS__
#define GPRINT(s) (std::string("[%f]\t") + s + std::string("\n")).c_str(), TPLUS(this->start_time)
#define FPRINTF(f, s, v, ...) { if(v) { std::fprintf(f, GPRINTF(s, __VA_ARGS__)); std::fflush(f); } }
#define FPRINT(f, s, v) { if(v) { std::fprintf(f, GPRINT(s)); std::fflush(f); } }
#define PRINTF(s, v, ...) { if(v) { std::printf(GPRINTF(s, __VA_ARGS__)); std::fflush(stdout); } }
#define PRINT(s, v) { if(v) { std::printf(GPRINT(s)); std::fflush(stdout); } }
#define WPRINTF(s, ...) { FPRINTF(this->work_log, s, IS(LOG_WORK), __VA_ARGS__) PRINTF(s, IS(VER_WORK), __VA_ARGS__) }
#define DPRINTF(s, ...) { FPRINTF(this->data_log, s, IS(LOG_DATA), __VA_ARGS__) PRINTF(s, IS(VER_DATA), __VA_ARGS__) }
#define WPRINT(s) { FPRINT(this->work_log, s, IS(LOG_WORK)) PRINT(s, IS(VER_WORK)) }
#define DPRINT(s) { FPRINT(this->data_log, s, IS(LOG_DATA)) PRINT(s, IS(VER_DATA)) }

// The rec_gen abstract class represents the general structure of
// the Rec-Gen algorithm; different child classes can implement
// the individual steps differently
class rec_gen
{
public:
	// Hypergraph member class for sibling detection
	class hypergraph
	{
	public:
		class edge {};
		virtual void insert_edge(edge e) {}
		virtual void erase_edge(edge e) {}
		virtual std::set<coupled_node*> extract_clique() {}
	};
protected:
	// Internal pedigree
	poisson_pedigree* ped;
	// Reconstruct the genetic material of top-level coupled node v (returns v)
	virtual coupled_node* collect_symbols(coupled_node* v) { return v; }
	// Perform statistical tests to detect siblinghood (returns hypergraph)
	virtual hypergraph* test_siblinghood() { return new hypergraph(); }
	// Assign parents to the top-level generation based on the siblinghood hypergraph
	virtual void assign_parents(hypergraph* G) {}
	// Logging info
	/// Log file names and pointers
	std::string work_path, data_path;
	std::FILE *work_log, *data_log;
	/// Time of reconstruction start
	std::chrono::high_resolution_clock::time_point start_time;
	/// Settings bitstring
	long long settings;
	// Initialize given all info
	void init(poisson_pedigree* ped, std::string work_log, std::string data_log, double sib, double rec, int d, long long settings);
	// Private members
	/// Search parameters
	double sib; /// Threshold of genetic overlap for sibling triples
	double rec; /// Proportion of genome that needs to be recovered for a node to be valid
	int d; /// The minimum desirable siblinghood clique size (Definition 4.2, d-richness)
public:
	// Constructors
	/// Given pedigree
	rec_gen(poisson_pedigree* ped);
	/// Given all
	rec_gen(poisson_pedigree* ped, std::string work_log, std::string data_log, double sib, double rec, int d, long long settings);
	// Initialize post-construction -- important if parameters like filenames changed since construction
	void init();
	// Access pedigree
	poisson_pedigree* get_pedigree();
	// Rebuild (returns reconstructed pedigree)
	virtual poisson_pedigree* apply_rec_gen();
};

#endif
