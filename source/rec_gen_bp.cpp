
/********************************************************************
* Defines belief propagation symbol-collection in REC-GEN
********************************************************************/

#include "rec_gen_bp.h"
#include "bp_message.h"

#include <cstring>
#include <cmath>

// The rec_gen_bp class implements belief-propagation symbol collection
// Override: reconstruct the genetic material of top-level coupled node v (returns v)
coupled_node* rec_gen_bp::collect_symbols(coupled_node *v)
{
	/// Initialize sets of genes
	v->all_des_genes = new std::unordered_set<gene>[this->ped->num_blocks()]();
	if (this->ped->all_genes == NULL) {
		this->ped->all_genes = new std::unordered_set<gene>[this->ped->num_blocks()]();
		for (auto ext : (*this->ped)[0]) {
			ext->all_des_genes = new std::unordered_set<gene>[this->ped->num_blocks()]();
			for (int b = 0; b < this->ped->num_blocks(); b++) {
				this->ped->all_genes[b].insert((*(*ext)[0])[b]);
				ext->all_des_genes[b].insert((*(*ext)[0])[b]);
			}
		}
	}
	/// Iterate over all blocks
	for (int b = 0; b < this->ped->num_blocks(); b++) {
		/// Find the set of all genes in subtree
		for (auto ext : *v)
			for (gene g : ext->couple()->all_des_genes[b])
				v->all_des_genes[b].insert(g);
		/// Initialize belief
		v->message(b, this->ped->all_genes->size(), std::pow(this->epsilon, v->num_ch()));
		/// Iterate over all pairs of genes
		for (gene g1 : v->all_des_genes[b])
			for (gene g2 : v->all_des_genes[b])
				if (g1 <= g2) {
					/// Set up DP
					long double num_missing_gene[v->num_ch() + 1][v->num_ch() + 1];
					std::memset(num_missing_gene, 0, sizeof num_missing_gene);
					num_missing_gene[0][0] = 1;
					/// DP over children
					int i = 1;
					for (individual_node* indiv : *v) {
						coupled_node* ch = indiv->couple();
						bp_message& message = *ch->message(b);
						for (int j = 0; j < v->num_ch(); j++) {
							long double p = message.get_marginal(g1) + (g1 != g2) * (message.get_marginal(g2) - message[bp_domain(g1, g2)]);
							num_missing_gene[i][j] += num_missing_gene[i - 1][j] * p;
							num_missing_gene[i][j + 1] += num_missing_gene[i - 1][j] * (1 - p);
						}
						i++;
					}
					/// Insert to distribution
					bp_message& message = *v->message(b);
					bp_domain domain = bp_domain(g1, g2);
					message.set(domain, 0);
					for (int i = 0; i <= v->num_ch(); i++)
						message.inc(domain, num_missing_gene[v->num_ch()][i] * std::pow(this->epsilon, i));
					WPRINTF("Marginal PDF (unnormalized) of genes %lld and %lld for couple %lld at block %d is %llf", g1, g2, v->get_id(), b, message[domain])
				}
		/// Normalize distribution
		v->message(b)->normalize();
		/// Pick maximum pair
		bp_domain genes = v->message(b)->extract_max();
		v->insert_gene(b, genes[0]), v->insert_gene(b, genes[1]);
		DPRINTF("For couple %lld at position %d found genes %lld and %lld (marginal: %llf)", v->get_id(), b, genes[0], genes[1], (*v->message(b))[genes])
	}
	/// Return argument node for chaining
	return v;
}

long double rec_gen_bp::set_epsilon(long double epsilon) { return this->epsilon = epsilon; }
