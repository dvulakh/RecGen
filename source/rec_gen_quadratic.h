
/********************************************************************
* An attempt to improve on the O((alpha^T)^3) runtime of
* rec_gen_basic.
********************************************************************/

#ifndef REC_GEN_QUADRATIC_H
#define REC_GEN_QUADRATIC_H

#include "rec_gen_basic.h"

// The rec_gen_basic class is a basic implementation of the Rec-Gen algorithm
// presented in the paper "Efficient Reconstruction of Stochastic Pedigrees"
class rec_gen_quadratic : public rec_gen_basic
{
protected:
	// Override: reconstruct the genetic material of top-level coupled node v (returns v)
	virtual coupled_node* collect_symbols(coupled_node* v);
	// Override: perform statistical tests to detect siblinghood (returns hypergraph)
	virtual hypergraph* test_siblinghood();
public:
	// Constructors
	/// Inherit
	using rec_gen_basic::rec_gen_basic;
};

#endif
