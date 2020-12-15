
/********************************************************************
* Implements recursive symbol-collection in REC-GEN
********************************************************************/

#include "rec_gen_recursive.h"
#include <algorithm>

// The rec_gen_recursive class implements recursive genome-finding
// Override: reconstruct the genetic material of top-level coupled node v (returns v)
coupled_node* rec_gen_recursive::collect_symbols(coupled_node* v)
{
	/// Iterate through all blocks
	for (int b = 0; b < this->ped->num_blocks(); b++) {
		/// Map genes to the vectors of bushiness values that occur
		std::unordered_map<gene, std::vector<int>> ch_block;
		for (individual_node* ch : *v)
			for (std::pair<gene, int> g : ch->couple()->get_des_genes(b)) {
				ch_block[g.first].push_back(g.second);
				DPRINTF("Found gene %lld at block %d of %lld for %lld, abundance %d", g.first, b, ch->couple()->get_id(), v->get_id(), g.second)
			}
		/// Sort children thresholds for each gene
		for (auto it = ch_block.begin(); it != ch_block.end(); it++)
			sort(it->second.begin(), it->second.end(), [&](int i, int j){ return i > j; });
		/// Determine the recursively satisfied bushiness for each gene
		/// Insert those genes that are above bush_th
		/// Also keep track of the best and second-best genes to use as guesses
		gene b1 = 0, b2 = 0;
		int th1 = 0, th2 = 0;
		for (auto it = ch_block.begin(); it != ch_block.end(); it++) {
			/// Get bushiness for this gene
			int th = 0;
			for (int i = 0; i < it->second.size(); i++)
				th = std::max(th, std::min(i + 1, it->second[i]));
			DPRINTF("At block %d of %lld, gene %lld has abundance %d", b, v->get_id(), it->first, th)
			/// Add if above threshold
			if (th >= this->bush_th)
				v->insert_des_gene(b, it->first, th);
			/// Consider for guesses
			if (th > th1) b2 = b1, th2 = th1, b1 = it->first, th1 = th;
			else if (th > th2) b2 = it->first, th2 = th;
		}
		/// Add guesses, inserting them also to the descendant genes list if necessary
		DPRINTF("For couple %lld at position %d found genes %lld and %lld (frequency: %d %d)", v->get_id(), b, b1, b2, th1, th2)
		v->insert_gene(b, b1)->insert_gene(b, b2);
		if (th1 < this->bush_th)
			v->insert_des_gene(b, b1, th1);
		if (th2 < this->bush_th)
			v->insert_des_gene(b, b2, th2);
	}
	return v;
}
/// Bushiness mutator
int rec_gen_recursive::set_bush_th(int bush_th)
{ return this->bush_th = bush_th; }
