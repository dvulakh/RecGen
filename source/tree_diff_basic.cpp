
/********************************************************************
* Implements a basic tree diff class that identifies original nodes
* to reconstructed nodes with similar children
********************************************************************/

#include "tree_diff_basic.h"

#include <algorithm>
#include <vector>

// Establish bijections of extant population based on id numbers
tree_diff* tree_diff_basic::biject_extant()
{
	/// Temporary ID map
	std::unordered_map<long long, individual_node*> orig_ids;
	for (coupled_node* coup : *this->orig)
		orig_ids[(*coup)[0]->get_id()] = (*coup)[0];
	/// Find same ids in reconstructed pedigree
	for (coupled_node* coup : *this->recon)
		this->biject(orig_ids[(*coup)[0]->get_id()]->couple(), coup);
	return this;
}

// Attempt to find same node in reconstructed tree based on children
// (assumes previous generation already bijected)
int tree_diff_basic::biject_parent(coupled_node* v)
{
	WPRINTF("Searching for image of %lldo in reconstructed pedigree", v->get_id())
	/// Best match: node pointer and number of correct children
	coupled_node* best = NULL;
	int num_match = 0;
	/// Get a map of possible parents to counts among children
	std::unordered_map<coupled_node*, int> par_count;
	/// Lambda for inserting parents to the map
	auto insert_to_pars = [&](coupled_node* par) {
		if (par) {
			par_count[par]++;
			DPRINTF("Adding node %lldr to parent counts: new count %d",
					par->get_id(), par_count[par])
		}
	};
	/// Loop through reconstructed children
	for (coupled_node* ch : *this->recon)
		if (this->re_to_or[ch] && v->is_child(this->re_to_or[ch])) {
			/// Try both parents
			insert_to_pars((*ch)[0]->parent());
			if ((*ch)[0] != (*ch)[1])
				insert_to_pars((*ch)[1]->parent());
		}
	/// Find the best candidate
	for (std::pair<coupled_node*, int> par : par_count)
		/// Make sure candidate satisfies minimum accuracy specs
		if (par.second > this->ch_acc * v->num_ch() &&
			par.second > this->ch_acc * par.first->num_ch() &&
			/// Make sure candidate is unclaimed
			(this->re_to_or.find(par.first) == this->re_to_or.end() || !this->re_to_or[par.first]) &&
			/// Make sure candidate is better than the current best
			par.second > num_match)
			/// Update best
			best = par.first, num_match = par.second;
	/// If no candidates, fail
	if (!best)
		return 0;
	/// If found a suitable candidate, insert
	this->biject(v, best);
	return num_match;
}

// Try to find a bijection between the topologies of the trees (return self)
tree_diff* tree_diff_basic::topology_biject()
{
	/// Set start time
	start_time = std::chrono::high_resolution_clock::now();
	/// Reset both pedigrees to extant layer
	this->orig->reset(), this->recon->reset();
	/// Establish bijections for the extant population
	this->biject_extant();
	/// Establish the bijections for ancestral generations
	while (!this->orig->done()) {
		/// Advance to next grade
		this->orig->next_grade();
		/// Add parents in this generation to vector and sort by number of children
		ADD_TO_BUCKET(nodes_total, this->orig->size());
		std::vector<std::pair<int, coupled_node*>> pars(this->orig->size());
		int i = 0;
		for (coupled_node* par : *this->orig)
			pars[i++] = { -par->num_ch(), par };
		std::sort(pars.begin(), pars.end());
		/// Find parent bijections
		for (std::pair<int, coupled_node*> par : pars) {
			/// Add number of edges
			ADD_TO_BUCKET(edges_total, par.second->num_ch());
			/// Try to find a match for the parent
			bool matched = this->biject_parent(par.second);
			/// Increment successful bijections
			ADD_TO_BUCKET(nodes_correct, matched);
			/// Count correct edges
			for (individual_node* ch : *par.second)
				if (matched && this->or_to_re[par.second]->is_child(this->or_to_re[ch->couple()]))
					ADD_TO_BUCKET(edges_correct, 1);
				else
					DPRINTF("Missing edge from (%lldo -> %lldr) to (%lldo -> %lldr)", par.second->get_id(), matched ? this->or_to_re[par.second]->get_id() : -1,
						ch->couple()->get_id(), this->or_to_re[ch->couple()] ? this->or_to_re[ch->couple()]->get_id() : -1)
		}
		/// Advance to next grade
		this->recon->next_grade();
	}
	/// Return self
	return this;
}

// Child accuracy threshold mutator
tree_diff_basic* tree_diff_basic::set_ch_acc(double ch_acc)
{ this->ch_acc = ch_acc; return this; }
