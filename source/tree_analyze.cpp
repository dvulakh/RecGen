
/********************************************************************
* Implements tools for analyzing the structure of poisson pedigrees:
* - Descendants of unique children of non-joint-LCA
* - Distribution of shared blocks in sibling and non-sibling triples
********************************************************************/

// TODO: Add inbreeding detection?

#include "tree_analyze.h"

// Perform preprocessing on tree
#define MERGE(x) merge(std::unordered_set<coupled_node*>(x))
preprocess::preprocess(poisson_pedigree* ped)
{
	// Claim pedigree
	this->ped = ped;
	// Associate vertices with their grade and descendants
	this->ped->reset();
	while (this->ped->cur_grade() < this->ped->num_grade()) {
		for (coupled_node *v : *this->ped) {
			/// Associate each vertex with the current grade
			this->grade_of[v] = this->ped->cur_grade();
			/// All nodes are in their own descendants pedigree
			this->des[v].insert(v);
			/// The extant generation are their own extant descendants
			if ((*v)[0] == (*v)[1]) {
				this->extant.push_back(v);
				this->ext[v].insert(v);
			}
			/// The remaining descendants of a node are the union of the descendants of its children
			for (individual_node* u : *v) {
				this->des[v].MERGE(this->des[u->couple()]);
				this->ext[v].MERGE(this->ext[u->couple()]);
			}
		}
		this->ped->next_grade();
	}
	// Associate vertices with their ancestors
	while(this->ped->cur_grade()) {
		this->ped->prev_grade();
		for (coupled_node* v : *this->ped) {
			/// All nodes are their own ancestors
			this->anc[v].insert(v);
			/// The remaining ancestors of a node are the ancestors of its parents
			for (int i = 0; i < 2; i++)
				this->anc[v].MERGE(this->anc[(*v)[i]->parent()]);
		}
	}
}

// Count for each generation the number of extant pair-vertex combinations
// such that the extant pair has a common ancestor that is a descendant of v
/// Returns a vector of pairs, where the first element at index i is the number
/// of such undesirable pairs for all v in generation i and the second is the
/// total number of pairs descended from unique children
std::vector<std::pair<long long, long long>> bad_joint_LCAs(preprocess* prep)
{
	/// Set up vector
	std::vector<std::pair<long long, long long>> bad_lca(prep->ped->num_grade(), {0,0});
	/// Iterate through all pairs of distinct extant nodes
	for (int x = 0; x < prep->extant.size(); x++)
		for (int y = x + 1; y < prep->extant.size(); y++) {
			/// Find their mutual ancestors
			std::unordered_set<coupled_node*> mut;
			for (coupled_node* v : prep->anc[prep->extant[x]])
				if (prep->anc[prep->extant[y]].find(v) != prep->anc[prep->extant[y]].end())
					mut.insert(v);
			/// Iterate over all pairs of mutual ancestors such that one is a descendant of another
			for (coupled_node* v : mut)
				for (coupled_node *u : mut)
					if (v != u && prep->des[v].find(u) != prep->des[v].end()) {
						/// If x and y are descended from unique children of v, this is a bad pair
						std::vector<individual_node*> chx, chy;
						std::unordered_set<coupled_node*> ancx = prep->anc[prep->extant[x]], ancy = prep->anc[prep->extant[y]];
						for (individual_node* ch : *v) {
							if (ancx.find(ch->couple()) != ancx.end())
								chx.push_back(ch);
							if (ancy.find(ch->couple()) != ancy.end())
								chy.push_back(ch);
						}
						if ((chx.size() > 1 || chy.size() > 1) || chx[0] != chy[0]) {
							bad_lca[prep->grade_of[v]].first++;
							break;
						}
					}
		}
	/// For each generation, find the number of pairs of elements descended from distinct children
	prep->ped->reset();
	while (!prep->ped->done()) {
		prep->ped->next_grade();
		/// For each vertex, get sum and sum of squares
		for (coupled_node* v : *prep->ped) {
			long long sig = 0, sig2 = 0;
			for (individual_node* ch : *v) {
				int nds = prep->ext[ch->couple()].size();
				sig += nds;
				sig2 += nds*nds;
			}
			bad_lca[prep->ped->cur_grade()].second += (sig*sig - sig2) / 2;
		}
	}
	return bad_lca;
}

// Report for each generation lists of the numbers of blocks shared by
// sibling and non-sibling triples
/// Returns a vector of triples of linked lists. Index 0 contains data for
/// unrelated triples; index 1 for tripling containing a sibling pair and an
/// unrelated third; index 2 for sibling triples
std::vector<std::list<int>*> block_share_stat(preprocess* prep)
{
	/// Set up vector
	std::vector<std::list<int>*> block_stat(prep->ped->num_grade(), NULL);
	for (std::list<int>*& e : block_stat)
		e = new std::list<int>[3]();
	/// For each generation, loop through all triples
	prep->ped->reset();
	while (!prep->ped->done()) {
		TRIPLE_IT(*prep->ped) {
			/// Determine how many pairs in the triple are siblings
			int nsib = 0, ns;
			coupled_node* ch[3] = { *u, *v, *w };
			for (int i = 0; i < 6; i++) {
				coupled_node* par = (*ch[i / 2])[i % 2]->parent();
				ns = 0;
				for (coupled_node* c : ch)
					ns += par->is_child(c);
				nsib = std::max(nsib, ns);
			}
			/// Insert shared block count
			block_stat[prep->ped->cur_grade()][std::max(0, nsib - 1)].push_back(shared_blocks(*u, *v, *w));
		}
		prep->ped->next_grade();
	}
	return block_stat;
}

// Report for each generation lists of the numbers of blocks shared by
// sibling triples
/// Returns a vector of linked lists: the distributions of shared blocks for
/// sibling triples in each generation
std::vector<std::list<int>> sib_block_share_stat(preprocess* prep)
{
	/// Set up vector
	std::vector<std::list<int>> block_stat(prep->ped->num_grade(), std::list<int>());
	/// For each generation loop through the triples of children of each couple
	prep->ped->reset();
	while (!prep->ped->done()) {
		prep->ped->next_grade();
		for (coupled_node* par : *prep->ped)
			TRIPLE_IT(*par)
				block_stat[prep->ped->cur_grade() - 1].push_back(shared_blocks((*u)->couple(), (*v)->couple(), (*w)->couple()));
	}
	return block_stat;
}
