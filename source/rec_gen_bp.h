
/********************************************************************
* Defines belief propagation symbol-collection in REC-GEN
********************************************************************/

#ifndef REC_GEN_BP_H
#define REC_GEN_BP_H

#include "rec_gen_quadratic.h"

#define MEM_PURGE_PAIRS (1 << 0)
#define MEM_PURGE_CHILD (1 << 1)

// The rec_gen_bp class implements belief-propagation symbol collection
class rec_gen_bp : public rec_gen_quadratic
{
protected:
	// Override: reconstruct the genetic material of top-level coupled node v (returns v)
	virtual coupled_node* collect_symbols(coupled_node* v);
	/// Compute one-time BP message helper
	bp_message& compute_message_at(coupled_node* v, int b);
	/// Probability assigned to event of finding a child with a gene not in its parents
	long double epsilon = 0.01;
	/// Types of strategies used to reduce memory footprint
	int memory_mode = 0;
public:
	/// Inherit constructor
	using rec_gen_quadratic::rec_gen_quadratic;
	/// Epsilon mutator
	long double set_epsilon(long double epsilon);
	/// Memory mode mutator
	int set_memory_mode(int memory_mode);
};

#endif
