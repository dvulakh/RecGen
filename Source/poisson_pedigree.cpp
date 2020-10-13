
/********************************************************************
* Implements the structures on which Rec-Gen operates: genes,
* individual nodes encapsulating the genetic information of a single
* person, coupled nodes representing mated pairs, and pedigree DAG
* of parent-child relationships
********************************************************************/

#include "poisson_pedigree.h"

#include <algorithm>
#include <random>
#include <vector>
#include <ctime>

INIT_ID(individual_node)
INIT_ID(coupled_node)

/************************** INDIVIDUALS ****************************/

// Initialize an individual node given all information
void individual_node::init(long long id, int genome_size, gene* genome, coupled_node* par, coupled_node* mate)
{
	id < 0 ? this->set_id() : this->set_id(id);
	this->genome_size = genome_size;
	this->genome = genome;
	this->par = par;
	this->mate = mate;
}

// Construct an individual node given the genome size and the ID
// For use during dump restoration
individual_node::individual_node(int genome_size, long long id)
{ init(id, genome_size, new gene[genome_size], NULL, NULL); }

// Construct an individual node given the genome size --
// initializes but does not fill genome
individual_node::individual_node(int genome_size)
{ init(-1, genome_size, new gene[genome_size], NULL, NULL); }

// Default constructor
individual_node::individual_node()
{ init(-1, 0, NULL, NULL, NULL); }

// Basic accessors
/// Return lvalue of block in index b
gene& individual_node::operator[](int b) { return this->genome[b]; }
/// Return the mate coupled node
coupled_node* individual_node::couple() { return this->mate; }
/// Get the parent couple
coupled_node* individual_node::parent() { return this->par; }

// Mate with another individual and return the couple
coupled_node* individual_node::mate_with(individual_node* other)
{
	this->mate = new coupled_node(this, other);
	return other->mate = this->mate;
}

// Assign a parent couple
coupled_node* individual_node::assign_par(coupled_node* par)
{ this->par = par; }

// Dump the individual information as a string
/// -i {id} -c {couple id} -p {parent id}
std::string individual_node::dump()
{
	std::string d = "-i " + std::to_string(this->get_id()) + " -c " +
		std::to_string(this->mate->get_id()) + " -p " +
		std::to_string(this->par->get_id()) + " -g " +
		std::to_string(this->genome_size);
	for (int i = 0; i < this->genome_size; i++)
		d += " " + std::to_string((*this)[i]);
	return d;
}

// Rebuild an individual from a dumped string
individual_node* individual_node::recover_dumped(std::string dump_out)
{
	// TODO
}

/**************************** COUPLES ******************************/

// Initialize an coupled node given all information
void coupled_node::init(long long id, std::pair<individual_node*, individual_node*> couple,
	std::unordered_set<individual_node*> children)
{
	id < 0 ? this->set_id() : this->set_id(id);
	this->couple = couple;
	this->children = children;
}

// Construct a coupled node given a pair to mate and the ID
// For use during dump restoration
coupled_node::coupled_node(individual_node* indiv1, individual_node* indiv2, long long id)
{ init(id, { indiv1, indiv2 }, std::unordered_set<individual_node*>()); }

// Construct a coupled node given a pair to mate
coupled_node::coupled_node(individual_node* indiv1, individual_node* indiv2)
{ init(-1, { indiv1, indiv2 }, std::unordered_set<individual_node*>()); }

// Construct a coupled node given an extant individual
coupled_node::coupled_node(individual_node* ext)
{ init(-1, { ext, ext }, std::unordered_set<individual_node*>()); }

// Access individuals by indexing
individual_node*& coupled_node::operator[](int index)
{
	if (index <= 0)
		return this->couple.first;
	return this->couple.second;
}

// Add a child to this couple's progeny
individual_node* coupled_node::add_child(individual_node* other)
{
	other->assign_par(this);
	this->children.insert(other);
	return other;
}

// Query if a node is a child
/// Query for an individual node
bool coupled_node::is_child(individual_node* other)
{
	/// Use the membership test of the map
	return this->children.find(other) != this->children.end();
}
/// Query for a coupled node
bool coupled_node::is_child(coupled_node* other)
{
	/// A coupled node is a child if either of its members is
	return this->is_child((*other)[0]) || this->is_child((*other)[1]);
}

// Iterating over a coupled_node iterates over its children
/// Internally, this is represented by iterating over the
/// children set
std::unordered_set<individual_node*>::iterator coupled_node::begin()
{ return this->children.begin(); }
std::unordered_set<individual_node*>::iterator coupled_node::end()
{ return this->children.end(); }

// Dump the couple information as a string
std::string coupled_node::dump()
{
	std::string d =  "-i " + std::to_string(this->get_id()) + " -m 2 " +
		   std::to_string((*this)[0]->get_id()) + " " +
		   std::to_string((*this)[1]->get_id()) + " -c " +
		   std::to_string(this->children.size());
	for (individual_node* ch : *this)
		d += " " + std::to_string(ch->get_id());
	return d;
}

// Rebuild a couple from a dumped string
coupled_node* coupled_node::recover_dumped(std::string dump_out)
{
	// TODO
}

/*********************** POISSON PEDIGREE **************************/

// Initialize a pedigree given all information
void poisson_pedigree::init(int genome_len, int tfr, int num_gen,
	int pop_sz, std::unordered_set<coupled_node*>* grades)
{
	this->genome_len = genome_len;
	this->tfr = tfr;
	this->num_gen = num_gen;
	this->cur_gen = -1;
	this->grades = grades;
	this->pop_sz = pop_sz;
	this->build();
}

// Build a poisson pedigree (4.1)
void poisson_pedigree::build()
{

	// Set up some temporary structures
	/// Temporarily store individual nodes in a vector as they're generated
	/// Assign a random value used for mating and use as sort key
	std::vector<std::pair<double, individual_node*>> mating_pool;
	/// The Poisson and uniform distributions used for generating fertility
	/// rate and mating individuals
	std::default_random_engine rng(time(NULL));
	std::poisson_distribution<int> poiss(this->tfr);
	std::uniform_real_distribution<double> unif(0, 1);

	// Generate the founder population
	this->cur_gen = this->num_gen - 1;
	/// If there is an odd member, ignore them, since they cannot mate
	for (int i = 0; i < this->pop_sz / 2 * 2; i++) {
		individual_node* indiv = new individual_node(this->genome_len);
		/// Give the ith founder gene i in all blocks
		for (int j = 0; j < this->genome_len; j++)
			(*indiv)[j] = i;
		/// Add to list with a random mating parameter
		mating_pool.push_back({ unif(rng), indiv });
	}

	// Mate the current generation, generate their children, and perform
	// symbol inheritance
	while (this->cur_gen) {
		/// Perform mating
		std::sort(mating_pool.begin(), mating_pool.end());
		for (int i = 0; i < mating_pool.size(); i += 2)
			this->add_to_current(mating_pool[i].second->mate_with(mating_pool[i + 1].second));
		mating_pool.clear();
		/// Generate children
		for (coupled_node* couple : *this) {
			int nch = poiss(rng);
			for (int i = 0; i < nch; i++) {
				individual_node* indiv = new individual_node(this->genome_len);
				/// Sample blocks randomly from parents
				for (int j = 0; j < this->genome_len; j++)
					(*indiv)[j] = (*(*couple)[unif(rng) < 0.5])[j];
				mating_pool.push_back({ unif(rng), couple->add_child(indiv) });
			}
		}
		/// Move down one grade
		this->cur_gen--;
	}

	// Couple the extant generation to itself
	for (auto pair_indiv : mating_pool)
		this->add_to_current(pair_indiv.second->mate_with(pair_indiv.second));

	// Set the founders as their own parents
	for (coupled_node* couple : this->grades[this->num_gen - 1])
		(*couple)[0]->assign_par(couple), (*couple)[1]->assign_par(couple);

}

// Construct given statistics, build a stochastic pedigree
poisson_pedigree::poisson_pedigree(int genome_len, int tfr, int num_gen, int pop_sz)
{ init(genome_len, tfr, num_gen, pop_sz, new std::unordered_set<coupled_node*>[num_gen]); }

// Statistic accessors
int poisson_pedigree::num_blocks() { return this->genome_len; }
int poisson_pedigree::num_child() { return this->tfr; }
int poisson_pedigree::num_grade() { return this->num_gen; }

// Adding and accessing coupled nodes in grades
/// Push an empty grade (returns self)
poisson_pedigree* poisson_pedigree::new_grade()
{
	/// Increment current generation counter and emplace a new set
	this->grades[++this->cur_gen] = std::unordered_set<coupled_node*>();
	return this;
}
/// Add to current grade (returns added node)
coupled_node* poisson_pedigree::add_to_current(coupled_node* couple)
{ this->grades[this->cur_gen].insert(couple); }
/// Access grades by indexing
std::unordered_set<coupled_node*>& poisson_pedigree::operator[](int grade)
{ return this->grades[grade]; }

// Iterating over a pedigree iterates over the last (highest)
// grade of coupled nodes
/// Internally, this is represented by iterating over the
/// current grade set
std::unordered_set<coupled_node*>::iterator poisson_pedigree::begin()
{ return this->grades[this->cur_gen].begin(); }
std::unordered_set<coupled_node*>::iterator poisson_pedigree::end()
{ return this->grades[this->cur_gen].end(); }

// Dump the pedigree information as a string
std::string poisson_pedigree::dump()
{
	// Get sets of all individuals and couples
	std::unordered_set<individual_node*> ind_set;
	std::unordered_set<coupled_node*> coup_set;
	/// Iterate over self at each generation
	for (this->cur_gen = 0; this->cur_gen < this->num_gen; this->cur_gen++)
		for (coupled_node* couple : *this)
			coup_set.insert(couple), ind_set.insert((*couple)[0]), ind_set.insert((*couple)[1]);
	// Dump
	std::string d;
	/// Dump individuals first
	for (individual_node* indiv : ind_set)
		d += "-i " + indiv->dump() + "\n";
	/// Then dump couples
	for (coupled_node* couple : coup_set)
		d += "-c " + couple->dump() + "\n";
	return d;
}

// Rebuild a pedigree from a dumped string
poisson_pedigree* poisson_pedigree::recover_dumped(std::string dump_out)
{
	// TODO
}
