
/********************************************************************
* Defines the structures on which Rec-Gen operates: genes, individual
* nodes encapsulating the genetic information of a single person,
* coupled nodes representing mated pairs, and pedigree DAG of parent-
* child relationships
********************************************************************/

#ifndef POISSON_PEDIGREE_H
#define POISSON_PEDIGREE_H

#include "flags.h"

#include <unordered_set>
#include <unordered_map>
#include <string>
#include <list>
#include <set>

struct individual_node;
struct coupled_node;
class poisson_pedigree;

struct bp_domain;

/********************** GENERAL DEFINITIONS ************************/

// Nodes are going to have IDs for the purpose of the dump/restore
// operation and for labelling isomorphisms
#define PRIVATE_ID_INFO(T) static long long ID_max; \
static std::unordered_map<long long, T*> ID_map; \
long long member_id; \
void set_id() { set_id(ID_max + 1); }
#define PUBLIC_ID_ACCESS(T) \
void set_id(long long id) { ID_map.erase(this->member_id); \
ID_map.insert({ this->member_id = id, this }); \
ID_max = ID_max < id ? id : ID_max; } \
static T* get_member_by_id(long long id) \
{ auto it = ID_map.find(id); return it == ID_map.end() ? NULL : it->second; } \
long long get_id() { return this->member_id; } \
static void clear_ids() { ID_map.clear(); ID_max = 0; }
#define INIT_ID(T) long long T::ID_max; std::unordered_map<long long, T*> T::ID_map;
#define NOT_COPYABLE(T) T(const T& other); T& operator=(const T&);

// Nodes and trees need to be dumpable and recoverable
#define DUMPABLE(T) std::string dump(); \
static T* recover_dumped(std::string dump_out, T*); \
static flag_reader frin;
#define INIT_DUMP(T) flag_reader T::frin;

// Iterating over all triples in L
#define TRIPLE_IT(L) for (auto u = (L).begin(); u != (L).end(); u++)\
for (auto v = std::next(u); v != (L).end(); v++)\
for (auto w = std::next(v); w != (L).end(); w++)

// Forward-declarations of classes to be imported
struct bp_message;

/************************** INDIVIDUALS ****************************/

// We represent a gene as a long long unsigned integer
// Permissible values: 0 - 18,446,744,073,709,551,615
typedef long long unsigned gene;

// Individual nodes encapsulate the genome of one person, the parent
// couple, and the mated coupled node
struct individual_node
{
private:
    // ID information
    PRIVATE_ID_INFO(individual_node)
    // Private members
    /// A genome is a sequence of blocks, each containing one gene
    /// The length of the gene array is determined by the pedigree
    gene* genome;
    int genome_size;
    /// The parent couple of this individual
    coupled_node* par;
    /// The coupled node this individual forms with their mate
    coupled_node* mate;
    // Private methods
    /// Initializer method chained from constructors
    void init(long long id, int genome_size, gene* genome, coupled_node* par, coupled_node* mate);
public:
    // No copying
    NOT_COPYABLE(individual_node)
    // Constructors
    /// Given the size of the genome and the id
    individual_node(int genome_size, long long id);
    /// Given the size of the genome
    individual_node(int genome_size);
    /// Default -- leaves genome NULL
    individual_node();
    // Destructor
    ~individual_node();
    individual_node* purge();
    // Accessors & mutators
    /// Index a modifiable lvalue of the gene in position b
    gene& operator[](int b);
    /// Get the genome size
    int num_blocks();
    /// Return the mate coupled node
    coupled_node* couple();
    /// Mate with another individual and return the couple
    coupled_node* mate_with(individual_node* other);
    /// Get the parent couple
    coupled_node* parent();
    /// Assign a parent couple
    coupled_node* assign_par(coupled_node* par);
    // Info dump
    DUMPABLE(individual_node)
    /// In addition to dumping full state,individual can dump just
    /// The id and genetic info
    std::string dump_genes();
    // ID Information
    PUBLIC_ID_ACCESS(individual_node)
};

/**************************** COUPLES ******************************/

// Coupled nodes encapsulate two individuals who have mated and
// their children
// As a special case, members of the unmated extant population are
// also considered coupled nodes
struct coupled_node
{
private:
    // ID information
    PRIVATE_ID_INFO(coupled_node)
    // Private members
    /// A coupled node contains one or two individual nodes
    /// Couples of one individual store two copies of that individual
    std::pair<individual_node*, individual_node*> couple;
    /// The children of a coupled node are internally stored in a
    /// hash table
    std::unordered_set<individual_node*> children;
    // Private methods
    /// Initializer method chained from constructors
    void init(long long id, std::pair<individual_node*, individual_node*> couple,
        std::unordered_set<individual_node*> children);
public:
    // No copying
    NOT_COPYABLE(coupled_node)
    // Constructors
    /// Given an id
    coupled_node(long long id);
    /// Given a pair to mate
    coupled_node(individual_node* indiv1, individual_node* indiv2);
    /// Given an extant individual
    coupled_node(individual_node* ext);
    /// Default
    coupled_node() {}
    // Destructor
    coupled_node* purge();
    // Access individuals by indexing
    individual_node*& operator[](int index);
    // Manipulate genes
    /// Query whether a member of the couple has gene g in block b
    bool has_gene(int b, gene g);
    /// Insert a gene to an unassigned couple member
    coupled_node* insert_gene(int b, gene g);
    // Get a parentless member of the couple
    individual_node* get_orphan();
    // Add a child to this couple's progeny
    individual_node* add_child(individual_node* other);
    // Query if a node is a child
    /// Query for an individual node
    bool is_child(individual_node* other);
    /// Query for a coupled node
    bool is_child(coupled_node* other);
    /// Query if two coupled nodes are siblings
    bool is_sib(coupled_node* other);
    /// Query for number of children
    int num_ch();
    // Remove child
    individual_node* erase_child(individual_node* ch);
    // Iterating over a coupled_node iterates over its children
    /// Internally, this is represented by iterating over the
    /// children set
    std::unordered_set<individual_node*>::iterator begin();
    std::unordered_set<individual_node*>::iterator end();
    // Get all extant descendants of a couple
    std::unordered_set<individual_node*> extant_desc(coupled_node* visitor=NULL);
    // Info dump
    DUMPABLE(coupled_node)
    // ID Information
    PUBLIC_ID_ACCESS(coupled_node)
// Extension for dfs pruning
private:
    /// Last vertex and block to visit this during dfs
    coupled_node *last_vis_vert;
// Extension for recursive symbol-collection
private:
    /// A linked list of genes this node "has" at each block
    /// First element is the gene itself; second is the minimum
    /// bushiness recursively satisfied
    std::list<std::pair<gene, int>>* rec_des_blocks;
    /// Initialize descendant gene array
    void init_des_blocks();
public:
    /// Get genes at block
    std::list<std::pair<gene, int>>& get_des_genes(int b);
    /// Insert a gene at block
    coupled_node* insert_des_gene(int b, gene g, int th);
// Extension for BP
private:
    /// Belief of gene distributions
    bp_message** belief;
    /// Size of genome
    int genome_len = -1;
    /// Last block accessed
    int last_block = -1;
public:
    /// Genes found in descendants pedigree
    std::unordered_set<gene>* all_des_genes;
    /// Belief accessor (initializes to current genes if NULL)
    bp_message*& message(int block, int domain_sz=0, long double nullval=0, int memory_mode=0);
// Extension for Parsimony
public:
    /// Set genes that participate in a minimum-error pair
    std::set<gene>* min_err;
};
// Count number of shared blocks in couple triple
int shared_blocks(coupled_node* u, coupled_node* v, coupled_node* w);

/*********************** POISSON PEDIGREE **************************/

// Pedigrees encapsulate the coupled nodes representing a population
// and organize other information, such as grades and growth rate
class poisson_pedigree
{
protected:
    // Private members
    /// Population information
    int genome_len; /// B in the paper
    int tfr; /// Expected number of children per couple; alpha
    int num_gen; /// T in the paper
    int pop_sz; /// N in the paper
    int cur_gen; /// The current grade number
    bool deterministic; /// Setting to true makes all fertilities
                        /// exactly alpha
    /// Grades of nodes are represented as unordered sets
    std::unordered_set<coupled_node*>* grades;
    // Private methods
    /// Initializer method chained from constructors
    void init(int genome_len, int tfr, int num_gen, int pop_sz, bool deterministic,
        std::unordered_set<coupled_node*>* grades);
public:
    // No copying
    NOT_COPYABLE(poisson_pedigree)
    // Constructors
    /// Given statistics, build a stochastic pedigree
    poisson_pedigree(int genome_len, int tfr, int num_gen, int pop_sz, bool deterministic);
    /// Default
    poisson_pedigree();
    /// Build pedigree
    poisson_pedigree* build();
    // Destructor
    poisson_pedigree* purge();
    // Statistic accessors
    int num_blocks();
    int num_child();
    int num_grade();
    int cur_grade();
    int size();
    // Adding and accessing coupled nodes in grades
    /// Reset the current grade pointer to zero (returns self)
    poisson_pedigree* reset();
    /// Return whether the current grade is the last one
    bool done();
    /// Push an empty grade (returns self)
    poisson_pedigree* new_grade();
    /// Move to the next/previous grade without pushing a new one (returns self)
    poisson_pedigree* next_grade();
    poisson_pedigree* prev_grade();
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
    // Info dump
    DUMPABLE(poisson_pedigree)
    /// In addition to dumping full info, a pedigree can dump just
    /// the extant population genetic data for REC-GEN input
    std::string dump_extant();
    // Generate pedigrees from shorthand
    static poisson_pedigree* parse_shorthand(std::string ped_string);
// Extension for BP
public:
    /// Set of all genes in extant population
    std::unordered_set<gene>* all_genes;
};

#endif
