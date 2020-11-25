
/********************************************************************
* Defines tools for analyzing the structure of poisson pedigrees:
* - Descendants of unique children of non-joint-LCA
********************************************************************/

// TODO: Add inbreeding detection?

#ifndef TREE_ANALYZE_H
#define TREE_ANALYZE_H

#include "poisson_pedigree.h"

// Structure that contains pedigree and preprocessing information
struct preprocess
{
	// The pedigree
	poisson_pedigree *ped;
	// Maps from each node to sets of ancestors and descendants
	std::unordered_map<coupled_node*, std::unordered_set<coupled_node*>> anc;
	std::unordered_map<coupled_node*, std::unordered_set<coupled_node*>> des;
	std::unordered_map<coupled_node*, std::unordered_set<coupled_node*>> ext;
	// Extant population vector
	std::vector<coupled_node*> extant;
	// Label each vertex with its generation
	std::unordered_map<coupled_node*, int> grade_of;
	// Constructor performs preprocessing
	preprocess(poisson_pedigree* ped);
};

// Count for each generation the number of extant pair-vertex combinations
// such that the extant pair has a common ancestor that is a descendant of v
/// Returns a vector of pairs, where the first element at index i is the number
/// of such undesirable pairs for all v in generation i and the second is the
/// total number of pairs descended from unique children
std::vector<std::pair<long long, long long>> bad_joint_LCAs(preprocess* prep);

#endif
