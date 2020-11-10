
/********************************************************************
* An attempt to improve on the O((alpha^T)^3) runtime of
* rec_gen_basic.
********************************************************************/

#include "rec_gen_quadratic.h"
#include <cstring>

#define NUM_BIT 32

/********************* QUADRATIC OVERRIDES *************************/

// Reconstruct the genetic material of top-level coupled node v (returns v)
/// IMPORTANT: All gene values must be no larger than extant population size!
coupled_node* rec_gen_quadratic::collect_symbols(coupled_node* par)
{
	// TODO: This step is making bad assumptions right now!

	/// Get a vector with all extant descendants by each direct child of par
	std::vector<std::unordered_set<individual_node*>> extant(par->num_ch());
	int i = 0; auto it = par->begin();
	for (; it != par->end(); i++, it++)
		extant[i] = (*it)->couple()->extant_desc();
	/// Populate for each block an array of which genes appear and with what frequency
	unsigned desc_have_gene[par->num_ch()][(*this->ped)[0].size() / NUM_BIT + 2];
	int num_block_appear[(*this->ped)[0].size() + 1];
	for (int b = 0; b < this->ped->num_blocks(); b++) {
		/// Reset the arrays
		std::memset(desc_have_gene, 0, sizeof(desc_have_gene));
		std::memset(num_block_appear, 0, sizeof(num_block_appear));
		/// For each direct child, collect the genes that appear among its extant descendants
		for (int ch = 0; ch < par->num_ch(); ch++) {
			for (individual_node *ext : extant[ch])
				desc_have_gene[ch][(*ext)[b] / NUM_BIT] |= 1 << (*ext)[b] % NUM_BIT;
			for (int g = 1; g <= (*this->ped)[0].size(); g++)
				num_block_appear[g] += (bool)(desc_have_gene[ch][g / NUM_BIT] & 1 << g % NUM_BIT);
		}
		/// Find the two most frequent genes in this block
		int g1 = 0, g2 = 0, c1 = 0, c2 = 0;
		for (int g = 1; g <= (*this->ped)[0].size(); g++)
			if (num_block_appear[g] > c1)
				g2 = g1, g1 = g, c2 = c1, c1 = num_block_appear[g];
			else if (num_block_appear[g] > c2)
				g2 = g, c2 = num_block_appear[g];
		/// Add those genes
		par->insert_gene(b, g1);
		par->insert_gene(b, g2);
	}
	/// Return the same couple
	return par;
}

// Perform statistical tests to detect siblinghood (returns hypergraph)
/// Currently considers all vertices, not just those with 99% rebuilt genomes
rec_gen::hypergraph* rec_gen_quadratic::test_siblinghood()
{
	/// Make a new graph
	rec_gen_quadratic::hypergraph_basic* G = new hypergraph_basic();
	/// Iterate over all pairs, inserting pairs that may form a triple into linked list
	DPRINT("Finding candidate pairs");
	std::list<std::pair<coupled_node*, coupled_node*>> sib_cand;
	for (auto it = ped->begin(); it != ped->end(); it++)
		for (auto jt = std::next(it); jt != ped->end(); jt++) {
			/// Count the number of shared blocks
			int shared_blocks = 0;
			for (int i = 0; i < this->ped->num_blocks(); i++)
				shared_blocks += (*jt)->has_gene(i, (*(**it)[0])[i]) || (*jt)->has_gene(i, (*(**it)[1])[i]);
			/// If the number of shared blocks is high enough, insert to candidates
			if (shared_blocks >= this->sib * this->ped->num_blocks())
				sib_cand.emplace_back(*it, *jt);
		}
	DPRINTF("Found %d candidate pairs; completing triples", sib_cand.size());
	/// For each pair, try to find a third element that completes the triple
	for (std::pair<coupled_node*, coupled_node*> pcc : sib_cand)
		for (coupled_node* coup : *this->ped)
			/// Make sure elements are distinct and triple has not yet been processed
			if (coup != pcc.first && coup != pcc.second && !G->query_edge({ pcc.first, pcc.second, coup })) {
				/// Count the number of shared blocks
				coupled_node *u = coup, *v = pcc.first, *w = pcc.second;
				int shared_blocks = 0;
				for (int i = 0; i < this->ped->num_blocks(); i++)
					shared_blocks += (v->has_gene(i, (*(*u)[0])[i]) && w->has_gene(i, (*(*u)[0])[i])) ||
									 (v->has_gene(i, (*(*u)[1])[i]) && w->has_gene(i, (*(*u)[1])[i]));
				/// If the number of shared blocks is high enough, insert a hyperedge
				if (shared_blocks >= this->sib * this->ped->num_blocks())
					G->insert_edge({ u, v, w });
				else DPRINTF("Failed triple completion: %lld %lld %lld", u->get_id(), v->get_id(), w->get_id());
			}
	/// Return the hypergraph
	return G;
}
