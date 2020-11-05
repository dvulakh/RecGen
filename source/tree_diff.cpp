
/********************************************************************
* Implements the general properties of a tree difference identifier
* for comparing the topologies and block contents of original and
* reconstructed pedigrees.
* Implementation details are filled in by child classes.
********************************************************************/

#include "tree_diff.h"

#include <cstring>

// Array resetting macro
#define RESET(arr) delete[] arr, arr = new int[this->orig->num_grade()], memset(arr, 0, sizeof(arr[0]) * this->orig->num_grade())

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
	// For each node that has an image in the reconstructed tree, compare their genomes
	this->orig->reset();
	while (!this->orig->done()) {
		this->orig->next_grade();
		ADD_TO_BUCKET(blocks_total, 2 * nodes_total_gen[this->orig->cur_grade()] * this->orig->num_blocks());
		for (coupled_node* coup : *this->orig) {
			/// Ignore extant nodes and NULL images
			auto it = this->or_to_re.find(coup);
			if (it == this->or_to_re.end() || !it->second)
				continue;
			auto bi = *it;
			if (bi.second && (*bi.first)[0] != (*bi.first)[1]) {
				int old_attempted = this->blocks_attempted, old_correct = this->blocks_correct;
				for (int b = 0; b < this->orig->num_blocks(); b++) {
					/// Add the number of non-zero blocks to the count of attempted blocks
					ADD_TO_BUCKET(blocks_attempted, (bool)(*(*bi.second)[0])[b] + (bool)(*(*bi.second)[1])[b]);
					/// Add the number of correct blocks
					ADD_TO_BUCKET(blocks_correct, (*(*bi.first)[0])[b] == (*(*bi.second)[0])[b] || (*(*bi.first)[0])[b] == (*(*bi.second)[1])[b]);
					/// If blocks are distinct, second block should match either
					ADD_TO_BUCKET(blocks_correct, (*(*bi.first)[1])[b] != (*(*bi.first)[0])[b] && ((*(*bi.first)[1])[b] == (*(*bi.second)[0])[b] || (*(*bi.first)[1])[b] == (*(*bi.second)[1])[b]));
					/// Otherwise, it should match both
					ADD_TO_BUCKET(blocks_correct, (*(*bi.first)[1])[b] == (*(*bi.first)[0])[b] && (*(*bi.first)[1])[b] == (*(*bi.second)[0])[b] && (*(*bi.first)[1])[b] == (*(*bi.second)[1])[b]);
				}
				DPRINTF("Comparing blocks in pair (%lldo -> %lldr): %d (%d%%) attempted; %d (%d%%/%d%%) correct", bi.first->get_id(), bi.second->get_id(),
					this->blocks_attempted - old_attempted, 50 * (this->blocks_attempted - old_attempted) / this->orig->num_blocks(),
					this->blocks_correct - old_correct, 50 * (this->blocks_correct - old_correct) / this->orig->num_blocks(),
					100 * (this->blocks_correct - old_correct) / std::max(1, this->blocks_attempted - old_attempted));
			}
		}
	}
	return this;
}

// Format a family of 7 ints into a statistics string
std::string tree_diff::stats_fmt(int node_t, int node_c, int edge_t, int edge_c, int block_t, int block_a, int block_c)
{
	node_t = std::max(node_t, 1), edge_t = std::max(edge_t, 1), block_t = std::max(block_t, 1), block_a = std::max(block_a, 1);
	return "Nodes correct:    " + std::to_string(node_c) + "/" + std::to_string(node_t) + " \t(" + std::to_string(100 * node_c / node_t) + "%)\n" +
		"Edges correct:    " + std::to_string(edge_c) + "/" + std::to_string(edge_t) + " \t(" + std::to_string(100 * edge_c / edge_t) + "%)\n" +
		"Blocks attempted: " + std::to_string(block_a) + "/" + std::to_string(block_t) + "\t(" + std::to_string(100 * block_a / block_t) + "%)\n" +
		"Blocks correct:   " + std::to_string(block_c) + "/" + std::to_string(block_t) + "\t(" + std::to_string(100 * block_c / block_t) + "%/" + std::to_string(100 * block_c / block_a) + "%)";
}

// Initialization and construction of tree_diff object
/// Initialize given all info
void tree_diff::init(poisson_pedigree* orig, poisson_pedigree* recon, std::string work_log, std::string data_log, long long settings)
{
	this->orig = orig;
	this->recon = recon;
	this->work_path = work_log;
	this->data_path = data_log;
	this->work_log = NULL;
	this->data_log = NULL;
	this->settings = settings;
	this->init();
}
/// Initialize based on current members
void tree_diff::init()
{
	/// Open files
	if (IS(LOG_WORK))
		delete this->work_log, this->work_log = std::fopen(this->work_path.c_str(), "w");
	if (IS(LOG_DATA))
		delete this->data_log, this->data_log = std::fopen(this->data_path.c_str(), "w");
	/// Create counter arrays based on pedigree height
	RESET(this->nodes_total_gen), RESET(this->nodes_correct_gen);
	RESET(this->edges_total_gen), RESET(this->edges_correct_gen);
	RESET(this->blocks_total_gen), RESET(this->blocks_attempted_gen), RESET(this->blocks_correct_gen);
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
