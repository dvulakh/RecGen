
/********************************************************************
* An attempt to improve on the O((alpha^T)^3) runtime of
* rec_gen_basic.
********************************************************************/

#ifndef REC_GEN_QUADRATIC_H
#define REC_GEN_QUADRATIC_H

#include "rec_gen_basic.h"

// The rec_gen_quadratic class improves the runtime of rec_gen_basic
// to quadratic in extant population size
class rec_gen_quadratic : public rec_gen_basic
{
protected:
    // Override: reconstruct the genetic material of top-level coupled node v (returns v)
    virtual coupled_node* collect_symbols(coupled_node* v);
    // Override: perform statistical tests to detect siblinghood (returns hypergraph)
    virtual hypergraph* test_siblinghood();
    // Whether to remove individuals from DFS consideration
    bool prune_dfs;
public:
    // Constructors
    /// Inherit
    using rec_gen_basic::rec_gen_basic;
    // Set DFS pruning
    rec_gen_quadratic* prune();
};

#endif
