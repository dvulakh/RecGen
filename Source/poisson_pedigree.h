
/********************************************************************
* Defines the structures on which Rec-Gen operates: genes, individual
* nodes encapsulating the genetic information of a single person,
* coupled nodes representing mated pairs, and pedigree DAG of parent-
* child relationships
********************************************************************/

#ifndef POISSON_PEDIGREE_H
#define POISSON_PEDIGREE_H

#include <unordered_set>
#include <string>

struct individual_node;
struct coupled_node;
class poisson_pedigree;

/************************** INDIVIDUALS ****************************/

// We represent a gene as a long long unsigned integer
// Permissible values: 0 - 18,446,744,073,709,551,615
typedef long long unsigned gene;

// Individual nodes encapsulate the genome of one person, the parent
// couple, and the mated coupled node
struct individual_node
{
private:
	// Private members
	/// A genome is a sequence of blocks, each containing one gene
	/// The length of the gene array is determined by the pedigree
	gene* genome;
	/// The parent couple of this individual
	coupled_node* par;
	/// The coupled node this individual forms with their mate
	coupled_node* mate;
	// Private methods
	/// Initializer method chained from constructors
	void init(gene* genome, coupled_node* par, coupled_node* mate);
public:
	// Constructors
	/// Given the size of the genome
	individual_node(int genome_size);
	/// Default -- leaves genome NULL
	individual_node();
	// Accessors & mutators
	/// Return a modifiable lvalue of the gene in position b
	gene& block(int b);
	/// Return the mate coupled node
	coupled_node* couple();
	/// Mate with another individual and return the couple
	coupled_node* mate_with(individual_node* other);
	/// Get the parent couple
	coupled_node* parent();
	/// Assign a parent couple
	coupled_node* assign_par(coupled_node* par);
};

/**************************** COUPLES ******************************/

// Coupled nodes encapsulate two individuals who have mated and
// their children
// As a special case, members of the unmated extant population are
// also considered coupled nodes
struct coupled_node
{
private:
	// Private members
	/// A coupled node contains one or two individual nodes
	/// Couples of one individual store two copies of that individual
	std::pair<individual_node*, individual_node*> couple;
	/// The children of a coupled node are internally stored in a
	/// hash table
	std::unordered_set<individual_node*> children;
	// Private methods
	/// Initializer method chained from constructors
	void init(std::pair<individual_node*, individual_node*> couple,
		std::unordered_set<individual_node*> children);
public:
	// Constructors
	/// Given a pair to mate
	coupled_node(individual_node* indiv1, individual_node* indiv2);
	/// Given an extant individual
	coupled_node(individual_node* ext);
	// Access individuals by indexing
	individual_node*& operator[](int index);
	// Add a child to this couple's progeny
	individual_node* add_child(individual_node* other);
	// Query if a node is a child
	/// Query for an individual node
	bool is_child(individual_node* other);
	/// Query for a coupled node
	bool is_child(coupled_node* other);
	// Iterating over a coupled_node iterates over its children
	/// Internally, this is represented by iterating over the
	/// children set
	std::unordered_set<individual_node*>::iterator begin();
	std::unordered_set<individual_node*>::iterator end();
};

/*********************** POISSON PEDIGREE **************************/

// Pedigrees encapsulate the coupled nodes representing a population
// and organize other information, such as grades and growth rate
class poisson_pedigree
{
private:
	// Private members
	/// Population information
	int genome_len; /// B in the paper
	int tfr; /// Expected number of children per couple; alpha
	int num_gen; /// T in the paper
	int cur_gen; /// The current grade number
	/// Grades of nodes are represented as unordered sets
	std::unordered_set<coupled_node*>* grades;
	// Private methods
	/// Initializer method chained from constructors
	void init(int genome_len, int tfr, int num_gen,
		std::unordered_set<coupled_node*>* grades);
	/// Build pedigree
	void build();
public:
	// Constructors
	/// Given statistics, build a stochastic pedigree
	poisson_pedigree(int genome_len, int tfr, int num_gen);
	// Statistic accessors
	int num_blocks();
	int num_child();
	int num_grade();
	// Adding and accessing coupled nodes in grades
	/// Push an empty grade (returns self)
	poisson_pedigree* new_grade();
	/// Add to current grade (returns added node)
	coupled_node* add_to_current(coupled_node* couple);
	/// Access grades by indexing
	std::unordered_set<coupled_node*>& operator[](int grade);
	// Iterating over a pedigree iterates over the last (highest)
	// grade of coupled nodes
	/// Internally, this is represented by iterating over the
	/// current grade set
	std::unordered_set<coupled_node*>::iterator begin();
	std::unordered_set<coupled_node*>::iterator end();
	// Dump the pedigree information as a string
	std::string dump();
};
// Rebuild a pedigree from a dumped string
poisson_pedigree* rebuild_dumped_pedigree(std::string dump);

#endif
