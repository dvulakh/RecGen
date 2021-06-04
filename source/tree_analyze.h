
/********************************************************************
* Defines tools for analyzing the structure of poisson pedigrees:
* - Descendants of unique children of non-joint-LCA
* - Distribution of shared blocks in sibling and non-sibling triples
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
    // Padding for printing
    int id_width, gene_width;
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

// Report for each generation lists of the numbers of blocks shared by
// sibling and non-sibling triples
/// Returns a vector of triples of linked lists. Index 0 contains data for
/// unrelated triples; index 1 for tripling containing a sibling pair and an
/// unrelated third; index 2 for sibling triples
std::vector<std::list<int>*> block_share_stat(preprocess* prep, int start);

// Report for each generation lists of the numbers of blocks shared by
// sibling triples
/// Returns a vector of linked lists: the distributions of shared blocks for
/// sibling triples in each generation
std::vector<std::list<int>> sib_block_share_stat(preprocess* prep);

// Display the induced pedigree of a given vertex, or the whole pedigree if
// the argument vertex is NULL
/// Prints out a depth-first traversal of the pedigree, replacing repeated
/// vertices with the message `[backedge]`
std::string print_sub_ped(preprocess* prep, coupled_node* v = NULL);

// Make tree-like pedigree
poisson_pedigree* tree_ped(int B, int T, int A);

#endif
