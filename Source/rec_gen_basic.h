
/********************************************************************
* Defines a basic implementation of the Rec-Gen algorithm.
********************************************************************/

#ifndef REC_GEN_BASIC_H
#define REC_GEN_BASIC_H

#include "rec_gen.h"

#include <initializer_list>
#include <map>
#include <set>

// The rec_gen_basic class is a basic implementation of the Rec-Gen algorithm
// presented in the paper "Efficient Reconstruction of Stochastic Pedigrees"
class rec_gen_basic : public rec_gen
{
public:
	// Basic hypergraph implementation
	class hypergraph_basic : public hypergraph
	{
	public:
		// Edge class
		class edge_basic
		{
		protected:
			// Couples contained in edge
			/// Internally sort by id so that a < b < c
			coupled_node *a, *b, *c;
		public:
			// Constructor
			/// Construct with initializer list
			edge_basic(std::initializer_list<coupled_node*> vert);
			// Iterator class
			/// Iterating over an edge iterates over contained vertices
			class iterator
			{
			protected:
				/// Index
				int i;
				/// Parent
				edge_basic& par;
			public:
				/// Construct
				iterator(int, edge_basic&);
				/// Advance
				iterator operator++();
				/// Compare
				bool operator==(const iterator&);
				bool operator!=(const iterator&);
				/// Dereference
				coupled_node* operator*();
				/// Referral
				coupled_node& operator->();
			};
			/// Iterators need access to parent privates
			friend class iterator;
			/// Begin and end iterators
			iterator begin();
			iterator end();
			// Comparison (for sets)
			bool operator<(const edge_basic&) const;
		};
	protected:
		// Graph information
		/// Vertex set -- set of all vertices and hyperdges that contain them
		std::map<coupled_node*, std::set<edge_basic>> vert;
		/// Adjacency map -- set of all edges and their multiplicities
		std::map<edge_basic, int> adj;
		/// Best clique found
		std::set<coupled_node*> best_clique;
		/// Clique currently under construction
		std::set<coupled_node*> clique;
		// Recursively build best clique
		void construct_best_clique_BB(std::map<coupled_node*, std::set<edge_basic>>::iterator it, int depth);
	public:
		// Constructor
		hypergraph_basic();
		// Inherited methods
		virtual void insert_edge(edge_basic e);
		virtual void erase_edge(edge_basic e);
		virtual std::set<coupled_node*> extract_clique();
	};
protected:
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
	rec_gen_basic(poisson_pedigree* ped, std::string work_log, std::string data_log, double sib, double rec, int d, long long settings) :
	rec_gen(ped, work_log, data_log, sib, rec, d, settings) {}
};

#endif
