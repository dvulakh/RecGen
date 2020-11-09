
/********************************************************************
* A basic implementation of the Rec-Gen algorithm.
********************************************************************/

#include "rec_gen_basic.h"

#define TRIPLE_IT(L) for (auto u = (L).begin(); u != (L).end(); u++)\
for (auto v = std::next(u); v != (L).end(); v++)\
for (auto w = std::next(v); w != (L).end(); w++)


/************************ BASIC REC-GEN ****************************/

// Reconstruct the genetic material of top-level coupled node v (returns v)
coupled_node* rec_gen_basic::collect_symbols(coupled_node* par)
{
	// TODO: This step is making bad assumptions right now!

	/// For all triples of immediate children of the parent couple
	TRIPLE_IT(*par)
		/// For all triples of distinct descendants of those childrenhs
		for (auto x : (*u)->couple()->extant_desc())
		for (auto y : (*v)->couple()->extant_desc())
		for (auto z : (*w)->couple()->extant_desc())
			if (x != y && y != z && z != x) {
				/// For all blocks
				DPRINTF("Identified extant triple: %lld %lld %lld", x->get_id(), y->get_id(), z->get_id());
				for (int i = 0; i < this->ped->num_blocks(); i++)
					/// If there is a shared value in this block, insert it
					if ((*x)[i] && (*x)[i] == (*y)[i] && (*y)[i] == (*z)[i] && !par->has_gene(i, (*x)[i])) {
						DPRINTF("Found gene for couple %lld (%lld %lld) at block %d: %lld", par->get_id(), (*par)[0]->get_id(), (*par)[1]->get_id(), i, (*x)[i]);
						par->insert_gene(i, (*x)[i]);
					}
			}
	/// Return same node
	return par;
}

// Perform statistical tests to detect siblinghood (returns hypergraph)
/// Currently considers all vertices, not just those with 99% rebuilt genomes
rec_gen::hypergraph* rec_gen_basic::test_siblinghood()
{
	/// Make a new graph
	rec_gen_basic::hypergraph_basic* G = new hypergraph_basic();
	/// Iterate over all triples
	TRIPLE_IT(*(this->ped)) {
		/// Count the number of shared blocks
		int shared_blocks = 0;
		for (int i = 0; i < this->ped->num_blocks(); i++)
			shared_blocks += ((*v)->has_gene(i, (*(**u)[0])[i]) && (*w)->has_gene(i, (*(**u)[0])[i])) ||
				((*v)->has_gene(i, (*(**u)[1])[i]) && (*w)->has_gene(i, (*(**u)[1])[i]));
		/// If the number of shared blocks is high enough, insert a hyperedge
		if (shared_blocks >= this->sib * this->ped->num_blocks())
			G->insert_edge({ *u, *v, *w });
	}
	/// Return the hypergraph
	return G;
}

// Assign parents to the top-level generation based on the siblinghood hypergraph
void rec_gen_basic::assign_parents(rec_gen::hypergraph* GARG)
{
	/// Convert to basic hypergraph
	hypergraph_basic* G = dynamic_cast<hypergraph_basic*>(GARG);
	/// Advance to next generation
	this->ped->new_grade();
	/// Repeatedly grab cliques and create parents for them
	std::set<coupled_node*> clique = G->extract_clique(d);
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
		clique = G->extract_clique(d);
		DPRINTF("Got a clique of size %d", clique.size());
	}
}

/********************** HYPERGRAPH STRUCTURE ***********************/

// HYPEREDGES

// Constructor -- use initializer list
rec_gen_basic::hypergraph_basic::edge_basic::edge_basic(std::initializer_list<coupled_node *> vert)
{
	/// Get values in the order that they appear
	auto it = vert.begin(); a = *(it++), b = (*it++), c = *(it++);
	/// Internally keep them sorted
	if (b->get_id() > c->get_id())
		std::swap(b, c);
	if (a->get_id() > b->get_id())
		std::swap(a, b);
	if (b->get_id() > c->get_id())
		std::swap(b, c);
}

// Iterator
/// Construct an iterator
rec_gen_basic::hypergraph_basic::edge_basic::iterator::iterator(int i, edge_basic& par) : par(par)
{ this->i = i; }
/// Advance an iterator
rec_gen_basic::hypergraph_basic::edge_basic::iterator rec_gen_basic::hypergraph_basic::edge_basic::iterator::operator++()
{
	iterator old = *this;
	this->i++;
	return old;
}
/// Compare two iterators
bool rec_gen_basic::hypergraph_basic::edge_basic::iterator::operator==(const iterator & ot)
{ return this->i == ot.i && !(this->par < ot.par) && !(ot.par < this->par); }
bool rec_gen_basic::hypergraph_basic::edge_basic::iterator::operator!=(const iterator & ot)
{ return !(*this == ot); }
/// Dereference an iterator
coupled_node* rec_gen_basic::hypergraph_basic::edge_basic::iterator::operator*()
{ return this->i > 0 ? this->i > 1 ? this->par.c : this->par.b : this->par.a; }
/// Referral from an iterator
coupled_node& rec_gen_basic::hypergraph_basic::edge_basic::iterator::operator->()
{ return ***this; }
/// Begin and end
rec_gen_basic::hypergraph_basic::edge_basic::iterator rec_gen_basic::hypergraph_basic::edge_basic::begin()
{ return iterator(0, *this); }
rec_gen_basic::hypergraph_basic::edge_basic::iterator rec_gen_basic::hypergraph_basic::edge_basic::end()
{ return iterator(3, *this); }

// Comparison (for sets)
bool rec_gen_basic::hypergraph_basic::edge_basic::operator<(const edge_basic& ot) const
{
	return this->a->get_id() == ot.a->get_id() ?
		this->b->get_id() == ot.b->get_id() ?
		this->c->get_id() < ot.c->get_id() :
		this->b->get_id() < ot.b->get_id() :
		this->a->get_id() < ot.a->get_id();
}

// GRAPH LOGIC

// Constructor -- create an empty hypergraph
rec_gen_basic::hypergraph_basic::hypergraph_basic()
{
	this->vert = std::map<coupled_node*, std::set<edge_basic>>();
	this->adj = std::map<edge_basic, int>();
}

// Insert an edge to the hypergraph
void rec_gen_basic::hypergraph_basic::insert_edge(edge_basic e)
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
void rec_gen_basic::hypergraph_basic::erase_edge(edge_basic e)
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
		if (v->get_orphan()->parent() != NULL) {
			for (edge_basic ve : this->vert[v])
				this->adj.erase(ve);
			this->vert.erase(v);
		}
}

// Extract a maximal hypergraph clique
/// Return whether vertex can be added to current clique
bool rec_gen_basic::hypergraph_basic::cliquable(coupled_node* vrt)
{
	for (auto u = clique.begin(); u != clique.end(); u++)
		for (auto v = std::next(u); v != clique.end(); v++)
			if (this->adj.find(edge_basic({ *u, *v, vrt })) == this->adj.end())
				return false;
	return true;
}
/// Find a clique of size d, if one exists
rec_gen_basic::hypergraph_basic* rec_gen_basic::hypergraph_basic::find_d_clique(
	std::map<coupled_node*, std::set<edge_basic>>::iterator it, int d)
{
	/// If there is already a clique of the necessary size, terminate
	if (this->clique.size() >= d)
		return this;
	/// Make sure the iterator is valid
	if (it == this->vert.end())
		return this;
	/// Try adding the current vertex
	if (this->cliquable(it->first)) {
		this->clique.insert(it->first);
		this->find_d_clique(std::next(it), d);
		/// Check whether resulting clique is satisfactory
		if (this->clique.size() >= d)
			return this;
		this->clique.erase(it->first);
	}
	/// If no d-size clique found, try adding the next vertex
	this->find_d_clique(std::next(it), d);
	return this;
}
/// Augment current clique so that it is maximal
rec_gen_basic::hypergraph_basic* rec_gen_basic::hypergraph_basic::augment_clique(
	std::map<coupled_node*, std::set<edge_basic>>::iterator it)
{
	/// Make sure the iterator is valid
	if (it == this->vert.end())
		return this;
	/// Try adding the current vertex
	if (this->cliquable(it->first))
		this->clique.insert(it->first);
	/// See if any other vertices can be added
	this->augment_clique(std::next(it));
	return this;
}
/// Do some setup and call the recursive implementation
std::set<coupled_node*> rec_gen_basic::hypergraph_basic::extract_clique(int d)
{
	/// Reset clique
	this->clique = std::set<coupled_node*>();
	/// Run recursion
	this->find_d_clique(this->vert.begin(), d)->augment_clique(this->vert.begin());
	/// Return clique
	return this->clique;
}
