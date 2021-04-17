
/********************************************************************
* Implements the structures on which Rec-Gen operates: genes,
* individual nodes encapsulating the genetic information of a single
* person, coupled nodes representing mated pairs, and pedigree DAG
* of parent-child relationships
********************************************************************/

#include "poisson_pedigree.h"
#include "rec_gen_bp.h"
#include "bp_message.h"

#include <algorithm>
#include <cstring>
#include <sstream>
#include <random>
#include <vector>
#include <ctime>

INIT_ID(individual_node)
INIT_ID(coupled_node)

INIT_DUMP(individual_node)
INIT_DUMP(coupled_node)
INIT_DUMP(poisson_pedigree)

/************************** INDIVIDUALS ****************************/

// Initialize an individual node given all information
void individual_node::init(long long id, int genome_size, gene* genome, coupled_node* par, coupled_node* mate)
{
	id < 0 ? this->set_id() : this->set_id(id);
	this->genome_size = genome_size;
	this->genome = genome;
	this->par = par;
	this->mate = mate;
	memset(this->genome, 0, sizeof(gene) * genome_size);
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

// Destructor
individual_node::~individual_node()
{
	/// Remove from children lists if present
	if (this->par != NULL)
		this->par->erase_child(this);
}
individual_node* individual_node::purge()
{
	/// Delete genome
	delete[] this->genome;
	return this;
}

// Basic accessors
/// Return lvalue of block in index b
gene& individual_node::operator[](int b) { return this->genome[b]; }
/// Return the genome size
int individual_node::num_blocks() { return this->genome_size; }
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
{ return this->par = par; }

// Dump the individual information as a string
/// -i {id} -c {couple id} -p {parent id} -g {genes}
std::string individual_node::dump()
{
	return "-i " + std::to_string(this->get_id()) + " -c " +
		std::to_string(this->mate ? this->mate->get_id() : 0) + " -p " +
		std::to_string(this->par ? this->par->get_id() : 0) + " " +
		this->dump_genes();
}
/// Dump only the genes: -i {id} -g {genes}
std::string individual_node::dump_genes()
{
	std::string d =  "-g " + std::to_string(this->genome_size);
	for (int i = 0; i < this->genome_size; i++)
		d += " " + std::to_string((*this)[i]);
	return d;
}

// Rebuild an individual from a dumped string
individual_node* individual_node::recover_dumped(std::string dump_out, individual_node* indiv)
{
	// Initialize the flag reader if not done already
	if (individual_node::frin.is_new()) {
		/// Read id
		frin.add_flag("id", 'i', 1, [&](std::vector<std::string> v, void* p) {
			individual_node* indiv = static_cast<individual_node*>(p);
			individual_node* orig = individual_node::get_member_by_id(std::stoll(v[0]));
			if (indiv == orig)
				return;
			individual_node::ID_map.erase(indiv->get_id());
			if (orig == NULL)
				indiv->set_id(std::stoll(v[0]));
			else {
				delete indiv;
				individual_node::frin.possess(orig);
			}
		});
		/// Read couple
		frin.add_flag("couple", 'c', 1, [&](std::vector<std::string> v, void* p) {
			static_cast<individual_node*>(p)->mate = coupled_node::get_member_by_id(std::stoll(v[0]));
		});
		/// Read parent
		frin.add_flag("parent", 'p', 1, [&](std::vector<std::string> v, void* p) {
			static_cast<individual_node*>(p)->assign_par(coupled_node::get_member_by_id(std::stoll(v[0])));
		});
		/// Read genome
		frin.add_flag("genome", 'g', -1, [&](std::vector<std::string> v, void* p) {
			individual_node* indiv = static_cast<individual_node*>(p);
			delete[] indiv->genome;
			indiv->genome_size = v.size();
			indiv->genome = new gene[indiv->genome_size];
			for (int i = 0; i < v.size(); i++)
				(*indiv)[i] = std::stoll(v[i]);
		});
	}
	// Possess the flag reader
	individual_node::frin.possess(indiv);
	individual_node::frin.read_flags(dump_out);
	return static_cast<individual_node*>(individual_node::frin.get_possessor());
}

/**************************** COUPLES ******************************/

// Initialize a coupled node given all information
void coupled_node::init(long long id, std::pair<individual_node*, individual_node*> couple,
	std::unordered_set<individual_node*> children)
{
	id < 0 ? this->set_id() : this->set_id(id);
	this->genome_len = -1;
	this->couple = couple;
	this->children = children;
	this->rec_des_blocks = NULL;
	this->all_des_genes = NULL;
	this->belief = NULL;
}

// Construct a coupled node given a pair to mate and the ID
// For use during dump restoration
coupled_node::coupled_node(long long id)
{ init(id, { NULL, NULL }, std::unordered_set<individual_node*>()); }

// Construct a coupled node given a pair to mate
coupled_node::coupled_node(individual_node* indiv1, individual_node* indiv2)
{ init(-1, { indiv1, indiv2 }, std::unordered_set<individual_node*>()); }

// Construct a coupled node given an extant individual
coupled_node::coupled_node(individual_node* ext)
{ init(-1, { ext, ext }, std::unordered_set<individual_node*>()); }

// Destructor for coupled node
coupled_node* coupled_node::purge()
{
	for (individual_node* ch : this->children)
		ch->assign_par(NULL);
	delete[] this->rec_des_blocks;
	delete[] this->all_des_genes;
	for (int b = 0; b <= genome_len; b++)
		delete this->belief[b];
	delete[] this->belief;
	return this;
}

// Access individuals by indexing
individual_node*& coupled_node::operator[](int index)
{
	if (index <= 0)
		return this->couple.first;
	return this->couple.second;
}

// Manipulate genes
/// Query whether a member of the couple has gene g in block b
bool coupled_node::has_gene(int b, gene g)
{ return g && ((*(*this)[0])[b] == g || (*(*this)[1])[b] == g); }
/// Insert a gene to an unassigned couple member
coupled_node* coupled_node::insert_gene(int b, gene g)
{ ((*(*this)[0])[b] ? *(*this)[1] : *(*this)[0])[b] = g; return this; }

// Get a parentless member of the couple
individual_node* coupled_node::get_orphan()
{ return (*this)[0]->parent() == NULL ? (*this)[0] : (*this)[1]; }

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
	/// NULL is not a child
	/// A coupled node is a child if either of its members is
	return other && (this->is_child((*other)[0]) || this->is_child((*other)[1]));
}
/// Query siblinghood
bool coupled_node::is_sib(coupled_node *other)
{
	/// Check if other is also a child of one of the parents
	return other && ((*this)[0]->parent()->is_child(other) || (*this)[1]->parent()->is_child(other));
}
/// Query for number of children
int coupled_node::num_ch() { return this->children.size(); }

// Remove a child from a couple's progeny
individual_node* coupled_node::erase_child(individual_node *ch)
{
	if (this->is_child(ch))
		this->children.erase(ch);
	return ch;
}

// Iterating over a coupled_node iterates over its children
/// Internally, this is represented by iterating over the
/// children set
std::unordered_set<individual_node*>::iterator coupled_node::begin()
{ return this->children.begin(); }
std::unordered_set<individual_node*>::iterator coupled_node::end()
{ return this->children.end(); }

// Get all extant descendants of a couple
std::unordered_set<individual_node*> coupled_node::extant_desc()
{
	/// If extant layer reached, return this
	if ((*this)[0] == (*this)[1])
		return std::unordered_set<individual_node*>({ (*this)[0] });
	/// Initialize a descendants set
	std::unordered_set<individual_node*> desc;
	/// For all children, add their descendants to this set
	for (individual_node* ch : (*this)) {
		auto ret = ch->couple()->extant_desc();
		desc.insert(ret.begin(), ret.end());
	}
	return desc;
}

// Dump the couple information as a string
/// -i {id} -m {member ids} -c {children ids}
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
coupled_node* coupled_node::recover_dumped(std::string dump_out, coupled_node* couple)
{
	// Initialize the flag reader if not done already
	if (coupled_node::frin.is_new()) {
		/// Read id
		frin.add_flag("id", 'i', 1, [&](std::vector<std::string> v, void* p) {
			coupled_node* couple = static_cast<coupled_node*>(p);
			coupled_node* orig = coupled_node::get_member_by_id(std::stoll(v[0]));
			if (couple == orig)
				return;
			coupled_node::ID_map.erase(couple->get_id());
			if (orig == NULL)
				couple->set_id(std::stoll(v[0]));
			else {
				delete couple;
				coupled_node::frin.possess(orig);
			}
		});
		/// Read members
		frin.add_flag("members", 'm', -1, [&](std::vector<std::string> v, void* p) {
			coupled_node* couple = static_cast<coupled_node*>(p);
			(*couple)[0] = individual_node::get_member_by_id(std::stoll(v[0]));
			(*couple)[1] = individual_node::get_member_by_id(std::stoll(v[1]));
		});
		/// Read children
		frin.add_flag("children", 'c', -1, [&](std::vector<std::string> v, void* p) {
			coupled_node* couple = static_cast<coupled_node*>(p);
			for (std::string s : v)
				couple->add_child(individual_node::get_member_by_id(std::stoll(s)));
		});
	}
	// Possess the flag reader
	coupled_node::frin.possess(couple);
	coupled_node::frin.read_flags(dump_out);
	return static_cast<coupled_node*>(coupled_node::frin.get_possessor());
}

// Extension for recursive symbol-collection
/// Initialize descendant gene array
void coupled_node::init_des_blocks()
{
	/// Create array
	this->rec_des_blocks = new std::list<std::pair<gene, int>>[(*this)[0]->num_blocks()]();
	/// Populate for extant node
	if ((*this)[0] == (*this)[1])
		for (int i = 0; i < (*this)[0]->num_blocks(); i++)
			this->rec_des_blocks[i].emplace_back((*(*this)[0])[i], INT32_MAX);
}
/// Get genes at block
std::list<std::pair<gene, int>>& coupled_node::get_des_genes(int b)
{
	/// Create blocks array on first query
	if (this->rec_des_blocks == NULL)
		this->init_des_blocks();
	/// Fetch required linked list
	return this->rec_des_blocks[b];
}
/// Insert a gene at block
coupled_node* coupled_node::insert_des_gene(int b, gene g, int th)
{
	/// Create blocks array on first query
	if (this->rec_des_blocks == NULL)
		this->init_des_blocks();
	/// Insert block
	this->rec_des_blocks[b].emplace_back(g, th);
	return this;
}

// Extension for belief-propagation
/// Return belief, initializing it to detected genes if null
bp_message*& coupled_node::message(int block, int domain_sz, long double nullval, int memory_mode)
{
	/// Initialize the belief if does not already exist
	if (!this->belief)
		/// If storing all of the messages, make an appropriate array
		if (!(memory_mode & MEM_PURGE_CHILD)) {
			this->genome_len = couple.first->num_blocks();
			this->belief = new bp_message*[couple.first->num_blocks()]();
			memset(this->belief, 0, sizeof(*this->belief) * this->genome_len);
		}
		/// Otherwise, make a size-1 array for the single gene
		else {
			this->belief = new bp_message*[1]();
			*this->belief = NULL;
		}
	/// Set the appropriate memory cell for mutation
	bp_message*& message_alias = memory_mode & MEM_PURGE_CHILD ? *this->belief : this->belief[block];
	/// Clear the message if it is stale
	if ((memory_mode & MEM_PURGE_CHILD) && this->last_block != block) {
		delete message_alias;
		message_alias = NULL;
	}
	/// If this message is computed, return it
	if (message_alias != NULL)
		return message_alias;
	/// If at an extant node, create a message with only those genes
	if (this->children.empty()) {
		message_alias = new bp_message(0, domain_sz);
		message_alias->inc(bp_domain((*(*this)[0])[block], (*(*this)[1])[block]), 1);
	}
	return message_alias;
}

// Count number of shared blocks in couple triple
int shared_blocks(coupled_node* u, coupled_node* v, coupled_node* w)
{
	int shr = 0;
	for (int i = 0; i < (*u)[0]->num_blocks(); i++)
		shr += (v->has_gene(i, (*(*u)[0])[i]) && w->has_gene(i, (*(*u)[0])[i])) ||
			(v->has_gene(i, (*(*u)[1])[i]) && w->has_gene(i, (*(*u)[1])[i]));
	return shr;
}

/*********************** POISSON PEDIGREE **************************/

// Initialize a pedigree given all information
void poisson_pedigree::init(int genome_len, int tfr, int num_gen, int pop_sz,
	bool deterministic, std::unordered_set<coupled_node*>* grades)
{
	this->genome_len = genome_len;
	this->tfr = tfr;
	this->num_gen = num_gen;
	this->cur_gen = -1;
	this->grades = grades ? grades : new std::unordered_set<coupled_node*>[num_gen];
	this->pop_sz = pop_sz;
	this->deterministic = deterministic;
	this->all_genes = NULL;
}

// Destructor: purge all pedigree members
poisson_pedigree* poisson_pedigree::purge()
{
	for (int grade = 0; grade < this->num_grade(); grade++)
		for (coupled_node* couple : this->grades[grade]) {
			for (individual_node* ch : *couple)
				delete ch->purge();
			if ((*couple)[0] != (*couple)[1])
				delete (*couple)[1]->purge();
			delete (*couple)[0]->purge();
			delete couple->purge();
		}
	delete[] this->grades;
	delete[] this->all_genes;
	return this;
}

// Build a poisson pedigree (4.1)
poisson_pedigree* poisson_pedigree::build()
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
	auto generate_fertility = [&](){ return this->deterministic ? this->tfr : poiss(rng); };

	// Generate the founder population
	this->cur_gen = this->num_gen - 1;
	/// If there is an odd member, ignore them, since they cannot mate
	for (int i = 1; i <= this->pop_sz / 2 * 2; i++) {
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
		/// Sort for mating
		std::sort(mating_pool.begin(), mating_pool.end());
		/// Delete last node of odd mating pool
		if (mating_pool.size() % 2) {
			delete mating_pool.back().second;
			mating_pool.pop_back();
		}
		for (int i = 0; i < mating_pool.size(); i += 2)
			this->add_to_current(mating_pool[i].second->mate_with(mating_pool[i + 1].second));
		mating_pool.clear();
		/// Generate children
		for (coupled_node* couple : *this) {
			int nch = generate_fertility();
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

	// Return self
	return this;

}

// Construct given statistics, build a stochastic pedigree
poisson_pedigree::poisson_pedigree(int genome_len, int tfr, int num_gen, int pop_sz, bool deterministic)
{ init(genome_len, tfr, num_gen, pop_sz, deterministic, new std::unordered_set<coupled_node*>[num_gen]); }

// Default constructor
poisson_pedigree::poisson_pedigree()
{ init(10, 3, 3, 10, 0, NULL); }

// Statistic accessors
int poisson_pedigree::num_blocks() { return this->genome_len; }
int poisson_pedigree::num_child() { return this->tfr; }
int poisson_pedigree::num_grade() { return this->num_gen; }
int poisson_pedigree::cur_grade() { return this->cur_gen; }
int poisson_pedigree::size() { return this->grades[this->cur_gen].size(); }

// Adding and accessing coupled nodes in grades
/// Reset the pedigree current grade pointer to 0 (returns self)
poisson_pedigree* poisson_pedigree::reset()
{
	this->cur_gen = 0;
	return this;
}
/// Return whether the current grade is the last one
bool poisson_pedigree::done() { return this->cur_gen == this->num_gen - 1; }
/// Push an empty grade (returns self)
poisson_pedigree* poisson_pedigree::new_grade()
{
	/// Increment current generation counter and emplace a new set
	this->grades[++this->cur_gen] = std::unordered_set<coupled_node*>();
	return this;
}
/// Move to the next/previous grade without pushing a new one (returns self)
poisson_pedigree* poisson_pedigree::next_grade()
{ this->cur_gen++; return this; }
poisson_pedigree* poisson_pedigree::prev_grade()
{ this->cur_gen--; return this; }
/// Add to current grade (returns added node)
coupled_node* poisson_pedigree::add_to_current(coupled_node* couple)
{ this->grades[this->cur_gen].insert(couple); return couple; }
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
	// Start with general info
	std::string d = "-B " + std::to_string(this->genome_len) +
			"\n-A " + std::to_string(this->tfr) +
			"\n-T " + std::to_string(this->num_gen) +
			"\n-N " + std::to_string(this->pop_sz) + "\n";
	// Get sets of all individuals and couples
	std::unordered_set<individual_node*> ind_set;
	std::unordered_set<coupled_node*> coup_set;
	/// Iterate over self at each generation
	for (this->cur_gen = 0; this->cur_gen < this->num_gen; this->cur_gen++)
		for (coupled_node* couple : *this)
			coup_set.insert(couple), ind_set.insert((*couple)[0]), ind_set.insert((*couple)[1]);
	// Dump nodes
	/// Prepare ids first
	for (individual_node* indiv : ind_set)
		d += "-i " + std::to_string(indiv->get_id()) + "\n";
	for (coupled_node* couple : coup_set)
		d += "-c " + std::to_string(couple->get_id()) + "\n";
	/// Then individuals
	for (individual_node* indiv : ind_set)
		d += "i " + indiv->dump() + "\n";
	/// Then dump couples
	for (coupled_node* couple : coup_set)
		d += "c " + couple->dump() + "\n";
	return d;
}

// Dump the extant population information as a string
std::string poisson_pedigree::dump_extant()
{
	// Dump the size of the extant population and the generation count
	std::string d = "-n " + std::to_string((*this)[0].size()) +
			"\n-T " + std::to_string(this->num_gen) +
			"\n-B " + std::to_string(this->genome_len) + "\n";
	// Dump the extant individual genetic data
	for (coupled_node* couple : (*this)[0])
		d += "i -i " + std::to_string((*couple)[0]->get_id()) + " " + (*couple)[0]->dump_genes() + "\n";
	return d;
}

// Rebuild a pedigree from a dumped string
poisson_pedigree* poisson_pedigree::recover_dumped(std::string dump_out, poisson_pedigree* ped)
{
	// Initialize the flag reader if not done already
	int extant_size = -1;
	if (poisson_pedigree::frin.is_new()) {
		/// Read block size
		frin.add_flag("blocks", 'B', 1, [&](std::vector<std::string> v, void* p) {
			static_cast<poisson_pedigree*>(p)->genome_len = std::stoi(v[0]);
		});
		/// Read alpha (TFR)
		frin.add_flag("alpha", 'A', 1, [&](std::vector<std::string> v, void* p) {
			static_cast<poisson_pedigree*>(p)->tfr = std::stoi(v[0]);
		});
		/// Read number of generations
		frin.add_flag("generations", 'T', 1, [&](std::vector<std::string> v, void* p) {
			poisson_pedigree* ped = static_cast<poisson_pedigree*>(p);
			ped->num_gen = std::stoi(v[0]);
			delete[] ped->grades;
			ped->grades = new std::unordered_set<coupled_node*>[ped->num_gen];
		});
		/// Read founder size
		frin.add_flag("founders", 'N', 1, [&](std::vector<std::string> v, void* p) {
			static_cast<poisson_pedigree*>(p)->pop_sz = std::stoi(v[0]);
		});
		/// Read extant size
		frin.add_flag("extant", 'n', 1, [&](std::vector<std::string> v, void* p) {
			extant_size = std::stoi(v[0]);
		});
		/// Read deterministic flag
		frin.add_flag("deterministic", 'd', 0, [&](std::vector<std::string> v, void* p) {
			static_cast<poisson_pedigree*>(p)->deterministic = true;
		});
		/// Make a new individual
		frin.add_flag("individual", 'i', 1, [&](std::vector<std::string> v, void* p) {
			new individual_node(1, std::stoll(v[0]));
		});
		/// Make a new couple
		frin.add_flag("couple", 'c', 1, [&](std::vector<std::string> v, void* p) {
			new coupled_node(std::stoll(v[0]));
		});
	}
	// Reset the identities of nodes
	individual_node::clear_ids();
	coupled_node::clear_ids();
	// Possess the flag reader
	poisson_pedigree::frin.possess(ped);
	// Prepare sets of nodes
	std::unordered_set<individual_node*> indivs;
	std::unordered_set<coupled_node*> coups;
	// Read lines
	std::istringstream sin(dump_out);
	std::string line;
	while(std::getline(sin, line)) {
		/// Ignore empty line
		if (line.empty())
			continue;
		/// Lines starting with '-' go directly to the pedigree flag reader
		else if (line[0] == '-')
			poisson_pedigree::frin.read_flags(line);
		/// Lines starting with 'i' go to the individual node dump restore
		else if (line[0] == 'i')
			indivs.insert(individual_node::recover_dumped(line, new individual_node()));
		/// Lines starting with 'c' go to the coupled node dump restore
		else if (line[0] == 'c')
			coups.insert(coupled_node::recover_dumped(line, new coupled_node()));
	}
	// If the extant generation parameter is set, insert the couples to the bottom
	if (extant_size >= 0) {
		/// Reset pedigree
		ped->reset();
		for (individual_node* indiv : indivs)
			ped->add_to_current(indiv->mate_with(indiv));
	}
	// Otherwise, find the extant population and rebuild the tree
	else {
		/// Set current generation to extant
		ped->reset();
		/// Extant is self-coupled
		for (coupled_node* couple : coups)
			if ((*couple)[0] == (*couple)[1])
				ped->add_to_current(couple);
		/// Add the ancestors
		while (ped->cur_gen < ped->num_gen - 1) {
			ped->new_grade();
			for (coupled_node* couple : (*ped)[ped->cur_gen - 1])
				for (int i = 0; i < 2; i++)
					if ((*couple)[i]->parent() && ped->grades[ped->cur_gen].find((*couple)[i]->parent()) == ped->end())
						ped->add_to_current((*couple)[i]->parent());
		}
	}
	return ped;
}
