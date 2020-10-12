
/********************************************************************
* Implements the structures on which Rec-Gen operates: genes,
* individual nodes encapsulating the genetic information of a single
* person, coupled nodes representing mated pairs, and pedigree DAG
* of parent-child relationships
********************************************************************/

#include "poisson_pedigree.h"

/************************** INDIVIDUALS ****************************/

// Initialize an individual node given all information
void individual_node::init(gene* genome, coupled_node* par, coupled_node* mate)
{
	this->genome = genome;
	this->par = par;
	this->mate = mate;
}

// Construct an individual node given the genome size --
// initializes but does not fill genome
individual_node::individual_node(int genome_size)
{ init(new gene[genome_size], NULL, NULL); }

// Default constructor
individual_node::individual_node()
{ init(NULL, NULL, NULL); }

// Basic accessors
/// Return lvalue of block in index b
gene& individual_node::block(int b) { return this->genome[b]; }
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

/**************************** COUPLES ******************************/

// Initialize an coupled node given all information
void coupled_node::init(std::pair<individual_node*, individual_node*> couple,
	std::unordered_set<individual_node*> children)
{
	this->couple = couple;
	this->children = children;
}

// Construct a coupled node given a pair to mate
coupled_node::coupled_node(individual_node* indiv1, individual_node* indiv2)
{ init({ indiv1, indiv2 }, std::unordered_set<individual_node*>()); }

// Construct a coupled node given an extant individual
coupled_node::coupled_node(individual_node* ext)
{ init({ ext, ext }, std::unordered_set<individual_node*>()); }

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

/*********************** POISSON PEDIGREE **************************/

// Initialize a pedigree given all information
void poisson_pedigree::init(int genome_len, int tfr, int num_gen,
	std::unordered_set<coupled_node*>* grades)
{
	this->genome_len = genome_len;
	this->tfr = tfr;
	this->num_gen = num_gen;
	this->cur_gen = -1;
	this->grades = grades;
	this->build();
}

// Build a poisson pedigree (4.1)
void poisson_pedigree::build()
{
	// TODO
}

// Construct given statistics, build a stochastic pedigree
poisson_pedigree::poisson_pedigree(int genome_len, int tfr, int num_gen)
{ init(genome_len, tfr, num_gen, new std::unordered_set<coupled_node*>[num_gen]); }

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
	// TODO
}
// Rebuild a pedigree from a dumped string
poisson_pedigree* rebuild_dumped_pedigree(std::string dump)
{
	// TODO
}
