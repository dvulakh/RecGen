
/********************************************************************
* Defines the general properties of a tree difference identifier for
* comparing the topologies and block contents of original and
* reconstructed pedigrees.
* Implementation details are filled in by child classes.
********************************************************************/

#ifndef TREE_DIFF_H
#define TREE_DIFF_H

#include "poisson_pedigree.h"
#include "logging.h"

#include <unordered_map>
#include <string>

// Help with counting statistics in both total and per-generation buckets
#define ADD_TO_BUCKET(arr, val) arr += val, arr##_gen[this->orig->cur_grade()] += val

// The tree_diff abstract class defines the operations that the checking
// program needs to support to assess rec_gen accuracy
class tree_diff
{
MAKE_LOGGABLE
protected:
	// The trees being compared
	poisson_pedigree *orig, *recon;
	// The bijection representing the topology
	std::unordered_map<coupled_node*, coupled_node*> or_to_re, re_to_or;
	// Add a matching of orig_vert to recon_vert in the bijection
	tree_diff* biject(coupled_node* orig_vert, coupled_node* recon_vert);
	// Initialize given all info
	void init(poisson_pedigree* orig, poisson_pedigree* recon, std::string work_log, std::string data_log, long long settings);
public:
	// Constructors
	/// Given two trees
	tree_diff(poisson_pedigree* orig, poisson_pedigree* recon);
	/// Given all
	tree_diff(poisson_pedigree* orig, poisson_pedigree* recon, std::string work_log, std::string data_log, long long settings);
	// Initialize post-construction -- important if parameters like filenames changed since construction
	void init();
	// Public information about diff results
	/// A full diff string
	std::string full_diff;
	/// Total number of nodes in orig tree, number correctly reconstructed
	int *nodes_total_gen = 0, *nodes_correct_gen = 0;
	int nodes_total, nodes_correct;
	/// Total number of edges in orig tree, number correctly reconstructed
	int *edges_total_gen = 0, *edges_correct_gen = 0;
	int edges_total, edges_correct;
	/// Total number of blocks in orig tree, number attempted/correct
	int *blocks_total_gen = 0, *blocks_attempted_gen = 0, *blocks_correct_gen = 0;
	int blocks_total, blocks_attempted, blocks_correct;
	/// Format a family of 7 ints into a statistics string
	static std::string stats_fmt(int node_t, int node_c, int edge_t, int edge_c, int block_t, int block_a, int block_c);
	// Try to find a bijection between the topologies of the trees (return self)
	virtual tree_diff* topology_biject() { return this; }
	// Check accuracy of assigned blocks once bijection is known (return self)
	virtual tree_diff* blocks_check();
};

#endif
