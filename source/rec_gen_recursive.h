
/********************************************************************
* Defines recursive symbol-collection in REC-GEN
********************************************************************/

#ifndef REC_GEN_RECURSIVE_H
#define REC_GEN_RECURSIVE_H

#include "rec_gen_quadratic.h"

// The rec_gen_recursive class implements recursive genome-finding
class rec_gen_recursive : public rec_gen_quadratic
{
protected:
	// Override: reconstruct the genetic material of top-level coupled node v (returns v)
	virtual coupled_node* collect_symbols(coupled_node* v);
	/// Minimum bushiness threshold for recursive symbol collection
	int bush_th;
public:
	/// Inherit constructor
	using rec_gen_quadratic::rec_gen_quadratic;
	/// Bushiness mutator
	int set_bush_th(int bush_th);
};

#endif
