
/********************************************************************
* Implements the structures on which Rec-Gen operates: genes,
* individual nodes encapsulating the genetic information of a single
* person, coupled nodes representing mated pairs, and pedigree DAG
* of parent-child relationships
********************************************************************/

#include "poisson_pedigree.h"

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
{ return this->par = par; }

// Dump the individual information as a string
/// -i {id} -c {couple id} -p {parent id} -g {genes}
std::string individual_node::dump()
{
	return "-i " + std::to_string(this->get_id()) + " -c " +
		std::to_string(this->mate->get_id()) + " -p " +
		std::to_string(this->par->get_id()) + " " +
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
	this->couple = couple;
	this->children = children;
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

	// Return self
	return this;

}

// Construct given statistics, build a stochastic pedigree
poisson_pedigree::poisson_pedigree(int genome_len, int tfr, int num_gen, int pop_sz)
{ init(genome_len, tfr, num_gen, pop_sz, new std::unordered_set<coupled_node*>[num_gen]); }

// Default constructor
poisson_pedigree::poisson_pedigree()
{ init(1, 1, 1, 1, NULL); }

// Statistic accessors
int poisson_pedigree::num_blocks() { return this->genome_len; }
int poisson_pedigree::num_child() { return this->tfr; }
int poisson_pedigree::num_grade() { return this->num_gen; }

// Adding and accessing coupled nodes in grades
/// Reset the pedigree current grade pointer to 0 (returns self)
poisson_pedigree* poisson_pedigree::reset()
{
	this->cur_gen = 0;
	return this;
}
/// Return whether the current grade is the last one
bool poisson_pedigree::done() { return this->cur_gen == this->num_gen; }
/// Push an empty grade (returns self)
poisson_pedigree* poisson_pedigree::new_grade()
{
	/// Increment current generation counter and emplace a new set
	this->grades[++this->cur_gen] = std::unordered_set<coupled_node*>();
	return this;
}
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
			"\n-T " + std::to_string(this->num_gen) + "\n";
	// Dump the extant individual genetic data
	for (coupled_node* couple : (*this)[0])
		d += "i " + (*couple)[0]->dump_genes() + "\n";
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
		/// Make a new individual
		frin.add_flag("individual", 'i', 1, [&](std::vector<std::string> v, void* p) {
			new individual_node(1, std::stoll(v[0]));
		});
		/// Make a new couple
		frin.add_flag("couple", 'c', 1, [&](std::vector<std::string> v, void* p) {
			new coupled_node(std::stoll(v[0]));
		});
	}
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
			ped->add_to_current(new coupled_node(indiv, indiv));
	}
	// Otherwise, find the founders and rebuild the tree
	else {
		/// Set current generation to founders
		ped->cur_gen = ped->num_gen - 1;
		/// Founders are their own parents
		for (coupled_node* couple : coups)
			if ((*couple)[0]->parent() == couple)
				ped->add_to_current(couple);
		/// Add the descendants
		ped->cur_gen--;
		while (ped->cur_gen >= 0) {
			for (coupled_node* couple : (*ped)[ped->cur_gen + 1])
				for (individual_node* ch : *couple)
					ped->add_to_current(ch->couple());
			ped->cur_gen--;
		}
	}
	return ped;
}
