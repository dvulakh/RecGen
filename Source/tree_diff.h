
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

// The tree_diff abstract class defines the operations that the checking
// program needs to support to assess rec_gen accuracy
class tree_diff
{
protected:
	MAKE_LOGGABLE
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
	int nodes_total, nodes_correct;
	/// Total number of edges in orig tree, number correctly reconstructed
	int edges_total, edges_correct;
	/// Total number of blocks in orig tree, number attempted/correct
	int blocks_total, blocks_attempted, blocks_correct;
	// Try to find a bijection between the topologies of the trees (return self)
	virtual tree_diff* topology_biject() { return this; }
	// Check accuracy of assigned blocks once bijection is known (return self)
	virtual tree_diff* blocks_check();
};

#endif
