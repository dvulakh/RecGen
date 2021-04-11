
/********************************************************************
* Defines belief propagation symbol-collection in REC-GEN
********************************************************************/

#ifndef REC_GEN_BP_H
#define REC_GEN_BP_H

#include "rec_gen_quadratic.h"

// The rec_gen_bp class implements belief-propagation symbol collection
class rec_gen_bp : public rec_gen_quadratic
{
protected:
	// Override: reconstruct the genetic material of top-level coupled node v (returns v)
	virtual coupled_node* collect_symbols(coupled_node* v);
	/// Probability assigned to event of finding a child with a gene not in its parents
	long double epsilon = 0.01;
	/// Whether to save space by purging pairs information after each round
	bool purge_pairs = false;
public:
	/// Inherit constructor
	using rec_gen_quadratic::rec_gen_quadratic;
	/// Epsilon mutator
	long double set_epsilon(long double epsilon);
	/// Memory mode mutator
	bool set_purge_pairs(bool purge_pairs);
};

#endif
