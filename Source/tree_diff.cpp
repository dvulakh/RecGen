
/********************************************************************
* Implements the general properties of a tree difference identifier
* for comparing the topologies and block contents of original and
* reconstructed pedigrees.
* Implementation details are filled in by child classes.
********************************************************************/

#include "tree_diff.h"

// Add a matching of orig_vert to recon_vert in the bijection
tree_diff* tree_diff::biject(coupled_node* orig_vert, coupled_node* recon_vert)
{
	DPRINTF("Established bijection (%lldo -> %lldr)", orig_vert->get_id(), recon_vert->get_id());
	this->or_to_re[orig_vert] = recon_vert;
	this->re_to_or[recon_vert] = orig_vert;
	return this;
}

// Count the number of attempted and correct blocks in the reconstructed tree
// (assuming a bijection is already present)
tree_diff* tree_diff::blocks_check()
{
	WPRINT(PRINT_HEADER("CHECKING BLOCKS"));
	// Set the total number of blocks based on the total number of nodes
	this->blocks_total = 2 * this->nodes_total * this->orig->num_blocks();
	// For each node that has an image in the reconstructed tree, compare their genomes
	for (std::pair<coupled_node*, coupled_node*> bi : this->or_to_re)
		/// Ignore extant nodes and NULL images
		if (bi.second && (*bi.first)[0] != (*bi.first)[1]) {
			int old_attempted = this->blocks_attempted, old_correct = this->blocks_correct;
			for (int b = 0; b < this->orig->num_blocks(); b++) {
				/// Add the number of non-zero blocks to the count of attempted blocks
				this->blocks_attempted += (bool)(*(*bi.second)[0])[b] + (bool)(*(*bi.second)[1])[b];
				/// Add the number of correct blocks
				this->blocks_correct += (*(*bi.first)[0])[b] == (*(*bi.second)[0])[b] || (*(*bi.first)[0])[b] == (*(*bi.second)[1])[b];
				/// If blocks are distinct, second block should match either
				this->blocks_correct += (*(*bi.first)[1])[b] != (*(*bi.first)[0])[b] && ((*(*bi.first)[1])[b] == (*(*bi.second)[0])[b] || (*(*bi.first)[1])[b] == (*(*bi.second)[1])[b]);
				/// Otherwise, it should match both
				this->blocks_correct += (*(*bi.first)[1])[b] == (*(*bi.first)[0])[b] && (*(*bi.first)[1])[b] == (*(*bi.second)[0])[b] && (*(*bi.first)[1])[b] == (*(*bi.second)[1])[b];
			}
			DPRINTF("Comparing blocks in pair (%lldo -> %lldr): %d (%d%%) attempted; %d (%d%%/%d%%) correct", bi.first->get_id(), bi.second->get_id(),
				this->blocks_attempted - old_attempted, 50 * (this->blocks_attempted - old_attempted) / this->orig->num_blocks(),
				this->blocks_correct - old_correct, 50 * (this->blocks_correct - old_correct) / this->orig->num_blocks(),
				100 * (this->blocks_correct - old_correct) / (this->blocks_attempted - old_attempted));
		}
	return this;
}

// Initialization and construction of tree_diff object
/// Initialize given all info
void tree_diff::init(poisson_pedigree* orig, poisson_pedigree* recon, std::string work_log, std::string data_log, long long settings)
{
	this->orig = orig;
	this->recon = recon;
	this->work_path = work_log;
	this->data_path = data_log;
	this->settings = settings;
	this->init();
}
/// Initialize based on current members
void tree_diff::init()
{
	/// Open files
	if (IS(LOG_WORK))
		this->work_log = std::fopen(this->work_path.c_str(), "w");
	if (IS(LOG_DATA))
		this->data_log = std::fopen(this->data_path.c_str(), "w");
	/// Clear counters
	this->nodes_total = 0;
	this->nodes_correct = 0;
	this->edges_total = 0;
	this->edges_correct = 0;
	this->blocks_total = 0;
	this->blocks_attempted = 0;
	this->blocks_correct = 0;

}
/// Construct given all info
tree_diff::tree_diff(poisson_pedigree* orig, poisson_pedigree* recon, std::string work_log, std::string data_log, long long settings)
{ init(orig, recon, work_log, data_log, settings); }
/// Construct given pedigree
tree_diff::tree_diff(poisson_pedigree* orig, poisson_pedigree* recon)
{ init(orig, recon, "tree-diff.log", "tree-dif.dat", ~0); }
