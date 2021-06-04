
/********************************************************************
* Defines parsimony for symbol-collection in REC-GEN
********************************************************************/

#ifndef REC_GEN_PARSIMONY_H
#define REC_GEN_PARSIMONY_H

#include "rec_gen_quadratic.h"

// The rec_gen_bp class implements belief-propagation symbol collection
class rec_gen_parsimony : public rec_gen_quadratic
{
protected:
    // Override: reconstruct the genetic material of top-level coupled node v (returns v)
    virtual coupled_node* collect_symbols(coupled_node* v);
public:
    /// Inherit constructor
    using rec_gen_quadratic::rec_gen_quadratic;
};

#endif
