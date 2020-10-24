
/********************************************************************
* A basic implementation of the Rec-Gen algorithm.
********************************************************************/

#include "rec_gen_basic.h"

#define TRIPLE_IT(L) for (auto u = (L).begin(); u != (L).end(); u++)\
for (auto v = std::next(u); v != (L).end(); v++)\
for (auto w = std::next(v); w != (L).end(); w++)

/********************** HYPERGRAPH STRUCTURE ***********************/

typedef std::set<coupled_node*> edge;

// Constructor -- create an empty hypergraph
rec_gen_basic::hypergraph_basic::hypergraph_basic()
{
	this->vert = std::map<coupled_node*, std::set<edge>>();
	this->adj = std::map<edge, int>();
}

// Insert an edge to the hypergraph
void rec_gen_basic::hypergraph_basic::insert_edge(edge e)
{
	/// Increment the number of occurrences of the edge by one
	this->adj[e]++;
	/// Maximum edge degree is 2 (per definition 3.11)
	this->adj[e] = std::max(this->adj[e], 2);
	/// Add all vertices to the vertex set
	for (coupled_node* v : e)
		this->vert[v].insert(e);
}

// Remove an edge from the hypergraph
void rec_gen_basic::hypergraph_basic::erase_edge(edge e)
{
	/// Decrement the number of occurrences of the edge by one
	this->adj[e]--;
	/// If no edges, erase map entry
	if (this->adj[e] <= 0) {
		this->adj.erase(e);
		/// Erase from all vertices
		for (coupled_node* v : e) {
			this->vert[v].erase(e);
			/// Erase vertex if no edges
			if (this->vert[v].empty())
				this->vert.erase(v);
		}
	}
	/// Also erase vertices if they have two assigned parents
	for (coupled_node* v : e)
		if (v->get_orphan()->parent() != NULL)
			this->vert.erase(v);
}

// Extract a maximal hypergraph clique
/// This method is naive recursive B&B implementation, but the author does not know
/// of polynomial-time algorithms for this task
void rec_gen_basic::hypergraph_basic::construct_best_clique_BB(
	std::map<coupled_node*, std::set<edge>>::iterator it, int depth)
{
	/// Check if a new best clique has been reached
	if (this->clique.size() > this->best_clique.size())
		this->best_clique = this->clique;
	/// Make sure the iterator is valid
	if (it == this->vert.end())
		return;
	/// Apply B&B -- if no possible clique reachable form this point is better than
	/// the one we have found so far, backtrack
	if (this->clique.size() + this->vert.size() - depth <= this->best_clique.size())
		return;
	/// Try not adding the the current vertex
	this->construct_best_clique_BB(std::next(it), depth + 1);
	/// See if it is possible to add the current vertex -- check all distinct pairs
	/// u, v already in the clique for hyperedges
	for (auto u = clique.begin(); u != clique.end(); u++)
		for (auto v = std::next(u); v != clique.end(); v++)
			if (this->adj.find(edge({ *u, *v, it->first })) == this->adj.end())
				return;
	/// Try adding the current vertex
	this->clique.insert(it->first);
	this->construct_best_clique_BB(std::next(it), depth + 1);
	/// Don't forget to remove before returning
	this->clique.erase(it->first);
}
/// Do some setup and call the recursive B&B implementation
std::set<coupled_node*> rec_gen_basic::hypergraph_basic::extract_clique()
{
	/// Reset cliques
	this->clique = this->best_clique = std::set<coupled_node*>();
	/// Run recursion
	this->construct_best_clique_BB(this->vert.begin(), 0);
	/// Return best clique
	return this->best_clique;
}

/************************ BASIC REC-GEN ****************************/

// Reconstruct the genetic material of top-level coupled node v (returns v)
coupled_node* rec_gen_basic::collect_symbols(coupled_node* par)
{
	// TODO: This step is making bad assumptions right now!

	/// For all triples of immediate children of the parent couple
	TRIPLE_IT(*par)
		/// For all triples of distinct descendants of those children
		for (auto x : (*u)->couple()->extant_desc())
		for (auto y : (*v)->couple()->extant_desc())
		for (auto z : (*w)->couple()->extant_desc())
			if (x != y && y != z && z != x) {
				/// For all blocks
				DPRINTF("Identified extant triple: %lld %lld %lld", x->get_id(), y->get_id(), z->get_id());
				for (int i = 0; i < this->ped->num_blocks(); i++)
					/// If there is a shared value in this block, insert it
					if ((*x)[i] && (*x)[i] == (*y)[i] && (*y)[i] == (*z)[i]) {
						DPRINTF("Found gene for couple %lld (%lld %lld) at block %d: %lld", par->get_id(), (*par)[0]->get_id(), (*par)[1]->get_id(), i, (*x)[i]);
						par->insert_gene(i, (*x)[i]);
					}
			}
	/// Return same node
	return par;
}

// Perform statistical tests to detect siblinghood (returns hypergraph)
rec_gen::hypergraph* rec_gen_basic::test_siblinghood()
{
	/// Make a new graph
	rec_gen_basic::hypergraph_basic* G = new hypergraph_basic();
	/// Iterate over all triples
	TRIPLE_IT(*(this->ped)) {
		/// Count the number of shared blocks
		int shared_blocks = 0;
		for (int i = 0; i < this->ped->num_blocks(); i++)
			shared_blocks += ((*v)->has_gene(i, (*(**u)[0])[i]) && (*w)->has_gene(i, (*(**u)[0])[i])) |
				((*v)->has_gene(i, (*(**u)[1])[i]) && (*w)->has_gene(i, (*(**u)[1])[i]));
		/// If the number of shared blocks is high enough, insert a hyperedge
		if (shared_blocks >= this->sib * this->ped->num_blocks())
			G->insert_edge({ *u, *v, *w });
	}
	/// Return the hypergraph
	return G;
}

// Assign parents to the top-level generation based on the siblinghood hypergraph
void rec_gen_basic::assign_parents(rec_gen::hypergraph* G)
{
	/// Advance to next generation
	this->ped->new_grade();
	/// Repeatedly grab cliques and create parents for them
	std::set<coupled_node*> clique = G->extract_clique();
	DPRINTF("Got a clique of size %d", clique.size());
	while (clique.size() >= d) {
		/// Create a new couple
		coupled_node* couple = (new individual_node(this->ped->num_blocks()))->mate_with(new individual_node(this->ped->num_blocks()));
		DPRINTF("Created new couple %lld (%lld, %lld)", couple->get_id(), (*couple)[0]->get_id(), (*couple)[1]->get_id());
		/// Add all of the clique elements to the couple's children
		for (coupled_node* ch : clique) {
			individual_node* orp = ch->get_orphan();
			DPRINTF("Attaching child %lld", orp->get_id());
			couple->add_child(orp);
		}
		/// Push the couple
		this->ped->add_to_current(couple);
		/// Remove one instance of each edge in the clique
		TRIPLE_IT(clique)
			G->erase_edge({ *u, *v, *w });
		/// Get the next clique
		clique = G->extract_clique();
		DPRINTF("Got a clique of size %d", clique.size());
	}
}
