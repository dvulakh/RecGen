
/********************************************************************
* Defines a basic tree diff class that identifies original nodes
* to reconstructed nodes with similar children
********************************************************************/

#ifndef TREE_DIFF_BASIC_H
#define TREE_DIFF_BASIC_H

#include "tree_diff.h"

// Default accuracy of child correspondences necessary to identify parent
#define DEFAULT_CH_ACC 0.49

// The tree_diff abstract class defines the operations that the checking
// program needs to support to assess rec_gen accuracy
class tree_diff_basic : public tree_diff
{
protected:
	// Percent of children that must be reconstructed for node to be considered reconstructed
	double ch_acc;
	// Establish bijections of extant population based on id numbers
	virtual tree_diff* biject_extant();
	// Attempt to find same node in reconstructed tree based on children
	// (assumes previous generation already bijected)
	// Returns number of matching children (0 if no match)
	virtual int biject_parent(coupled_node* v);
public:
	// Constructors
	/// Given two trees
	tree_diff_basic(poisson_pedigree* orig, poisson_pedigree* recon) :
	tree_diff(orig, recon) { this->ch_acc = DEFAULT_CH_ACC; }
	/// Given all
	tree_diff_basic(poisson_pedigree* orig, poisson_pedigree* recon, double ch_acc, std::string work_log, std::string data_log, long long settings) :
	tree_diff(orig, recon, work_log, data_log, settings) { this->ch_acc = ch_acc; }
	// Try to find a bijection between the topologies of the trees (return self)
	virtual tree_diff* topology_biject();
	// Mutator
	tree_diff_basic* set_ch_acc(double ch_acc);
};

#endif
