
/********************************************************************
* Defines the general structure of the Rec-Gen algorithm.
* Implementation details are filled in by child classes.
********************************************************************/

#ifndef REC_GEN_H
#define REC_GEN_H

#include "poisson_pedigree.h"
#include "logging.h"

#include <set>

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
	MAKE_LOGGABLE
	// Internal pedigree
	poisson_pedigree* ped;
	// Reconstruct the genetic material of top-level coupled node v (returns v)
	virtual coupled_node* collect_symbols(coupled_node* v) { return v; }
	// Perform statistical tests to detect siblinghood (returns hypergraph)
	virtual hypergraph* test_siblinghood() { return new hypergraph(); }
	// Assign parents to the top-level generation based on the siblinghood hypergraph
	virtual void assign_parents(hypergraph* G) {}
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
