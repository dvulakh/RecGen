
/********************************************************************
* Defines a basic implementation of the Rec-Gen algorithm.
********************************************************************/

#ifndef REC_GEN_BASIC_H
#define REC_GEN_BASIC_H

#include "rec_gen.h"

#include <map>

// The rec_gen_basic class is a basic implementation of the Rec-Gen algorithm
// presented in the paper "Efficient Reconstruction of Stochastic Pedigrees"
class rec_gen_basic : rec_gen
{
protected:
	// Basic hypergraph implementation
	class hypergraph_basic : hypergraph
	{
	protected:
		std::map<std::unordered_set<individual_node*>, int> adj;
	public:
		virtual void insert_edge(std::unordered_set<individual_node*> e);
		virtual void erase_edge(std::unordered_set<individual_node*> e);
		virtual std::unordered_set<individual_node*> extract_clique();
	};
	// Reconstruct the genetic material of top-level coupled node v (returns v)
	virtual coupled_node* collect_symbols(coupled_node* v);
	// Perform statistical tests to detect siblinghood (returns hypergraph)
	virtual hypergraph* test_siblinghood();
	// Assign parents to the top-level generation based on the siblinghood hypergraph
	virtual void assign_parents(hypergraph* G);
public:
	// Constructors
	/// Given pedigree
	rec_gen_basic(poisson_pedigree* ped) : rec_gen(ped) {}
	/// Given all
	rec_gen_basic(poisson_pedigree* ped, std::string work_log, std::string data_log, long long settings) :
	rec_gen(ped, work_log, data_log, settings) {}
};

#endif
