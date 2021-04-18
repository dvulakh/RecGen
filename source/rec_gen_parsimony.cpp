
/********************************************************************
* Implements parsimony for symbol-collection in REC-GEN
********************************************************************/

#include "rec_gen_parsimony.h"
#include "bp_message.h"

// Override: reconstruct the genetic material of top-level coupled node v (returns v)
coupled_node* rec_gen_parsimony::collect_symbols(coupled_node* v)
{
	/// Initialize the best-pairs map of v
	WPRINTF("Initializing min error sets for couple %lld", v->get_id())
	v->min_err = new std::set<gene>[this->ped->num_blocks()]();
	/// If children have uninitialized min_err, initialize to just genes
	for (individual_node* indiv : *v) {
		coupled_node *ch = indiv->couple();
		if (ch->min_err == NULL) {
			WPRINTF("Initializing min error sets for couple %lld", ch->get_id())
			ch->min_err = new std::set<gene>[this->ped->num_blocks()]();
			for (int b = 0; b < this->ped->num_blocks(); b++)
				ch->min_err[b].insert((*indiv)[b]);
		}
	}
	/// Process each block
	for (int b = 0; b < this->ped->num_blocks(); b++) {
		/// Get a set of all genes worth considering
		std::unordered_set<gene> des_genes;
		for (individual_node* indiv : *v)
			for (gene g : indiv->couple()->min_err[b])
				des_genes.insert(g);
		/// Iterate over all pairs, maintaining best current cost
		int best_cost = v->num_ch();
		std::set<bp_domain> min_pairs;
		for (gene g1 : des_genes)
			for (gene g2 : des_genes) {
				int cost = 0;
				for (individual_node* indiv : *v)
					cost += indiv->couple()->min_err[b].find(g1) == indiv->couple()->min_err[b].end() &&
						indiv->couple()->min_err[b].find(g2) == indiv->couple()->min_err[b].end();
				if (cost < best_cost)
					best_cost = cost, min_pairs.clear();
				if (cost == best_cost)
					min_pairs.insert(bp_domain(g1, g2));
				WPRINTF("Cost of genes %lld and %lld for couple %lld at block %d is %d", g1, g2, v->get_id(), b, cost)
			}
		/// Select an arbitrary pair of genes from the best
		v->insert_gene(b, (*min_pairs.begin())[0]);
		v->insert_gene(b, (*min_pairs.begin())[1]);
		DPRINTF("For couple %lld at position %d found genes %lld and %lld (cost: %d)", v->get_id(), b, (*(*v)[0])[b], (*(*v)[1])[b], best_cost)
		/// Collect all genes into min_err
		for (bp_domain d : min_pairs)
			v->min_err[b].insert(d[0]), v->min_err[b].insert(d[1]);
	}
	return v;
}
