
/********************************************************************
* Implements tools for analyzing the structure of poisson pedigrees:
* - Descendants of unique children of non-joint-LCA
* - Distribution of shared blocks in sibling and non-sibling triples
********************************************************************/

// TODO: Add inbreeding detection?

#include "tree_analyze.h"

#include <algorithm>

// Printing macros
#define UP_WIDTH(val, width) width = std::max(width, (int)std::to_string(val).length())
#define PAD(val, width) (std::string(width - std::to_string(val).length(), ' ') + std::to_string(val))
#define PAD_GENE(g) PAD(g, prep->gene_width)
#define PAD_ID(id) PAD(id, prep->id_width)

// Perform preprocessing on tree
#define MERGE(x) merge(std::unordered_set<coupled_node*>(x))
preprocess::preprocess(poisson_pedigree* ped)
{
	// Initial values
	this->id_width = this->gene_width = 0;
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
			/// Get widths for padding
			UP_WIDTH(v->get_id(), this->id_width);
			for (int b = 0; b < this->ped->num_blocks(); b++)
				UP_WIDTH((*(*v)[0])[b], this->gene_width), UP_WIDTH((*(*v)[1])[b], this->gene_width);
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
std::vector<std::list<int>*> block_share_stat(preprocess* prep, int start)
{
	/// Set up vector
	std::vector<std::list<int>*> block_stat(prep->ped->num_grade(), NULL);
	for (std::list<int>*& e : block_stat)
		e = new std::list<int>[3]();
	/// For each generation, loop through all triples
	prep->ped->reset();
	while (!prep->ped->done()) {
		if (prep->ped->cur_grade() >= start)
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

// Display the induced pedigree of a given vertex, or the whole pedigree if
// the argument vertex is NULL
/// Prints out a depth-first traversal of the pedigree, replacing repeated
/// vertices with the message `[backedge]`
std::string print_sub_ped(preprocess* prep, coupled_node* v)
{
	/// Recursion stack
	std::list<std::pair<coupled_node*, int>> stack;
	/// Answer string
	std::string ans;
	/// Visited couples
	std::unordered_set<coupled_node*> vis;
	/// Push v if it exists, otherwise the entire root population
	if (v) stack.push_back({ v, 0 });
	else
		for(coupled_node* couple : (*prep->ped)[prep->ped->num_grade() - 1])
			stack.push_back({ couple, 0 });
	/// For each couple on the queue, if it is not visited, expand and print it
	while (!stack.empty()) {
		/// Pop and print the entry
		auto front = stack.back(); stack.pop_back();
		std::string id_record = std::string(front.second, '>') + " " + PAD_ID(front.first->get_id()) +
			std::string(prep->ped->num_grade() - front.second, ' ');
		ans += id_record;
		/// If already visited, don't recurse
		if (vis.find(front.first) != vis.end()) {
			ans += "[backedge]\n";
			continue;
		}
		vis.insert(front.first);
		/// Print the genes
		ans += "|";
		for (int b = 0; b < prep->ped->num_blocks(); b++)
			ans += PAD_GENE((*(*front.first)[0])[b]) + " " + PAD_GENE((*(*front.first)[1])[b]) + " | ";
		ans[ans.length() - 1] = '\n';
		/// Push children
		for (individual_node* ch : *front.first)
			stack.emplace_back(ch->couple(), front.second + 1);
	}
	return ans;
}

// Generate a tree-pedigree
void tree_node(int, int, int, coupled_node*, gene&, poisson_pedigree*);
poisson_pedigree* tree_ped(int B, int T, int A)
{
	poisson_pedigree* ped = new poisson_pedigree(B, A, T, 0, 1);
	coupled_node* root = (new individual_node(B))->mate_with(new individual_node(B));
	(*ped)[T - 1].insert(root);
	for (int b = 0; b < B; b++)
		(*(*root)[0])[b] = 1, (*(*root)[1])[b] = 2;
	gene G = 3;
	for (int i = 0; i < A; i++)
		tree_node(B, T - 2, A, root, G, ped);
	return ped;
}
void tree_node(int B, int T, int A, coupled_node* par, gene& G, poisson_pedigree* ped)
{
	if (T) {
		coupled_node* couple = (new individual_node(B))->mate_with(par->add_child(new individual_node(B)));
		//coupled_node* couple = new coupled_node(new individual_node(B), par->add_child(new individual_node(B)));
		(*ped)[T].insert(couple);
		for (int b = 0; b < B; b++)
			(*(*couple)[0])[b] = G, (*(*couple)[1])[b] = (*(*par)[rand() > RAND_MAX / 2])[b];
		G++;
		for (int i = 0; i < A; i++)
			tree_node(B, T - 1, A, couple, G, ped);
	}
	else {
		individual_node* ext = par->add_child(new individual_node(B));
		for (int b = 0; b < B; b++)
			(*ext)[b] = (*(*par)[rand() > RAND_MAX / 2])[b];
		(*ped)[T].insert(ext->mate_with(ext));
	}
}
