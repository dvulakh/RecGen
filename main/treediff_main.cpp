
/********************************************************************
* Reads two poisson pedigrees (an original and a reconstructed
* version) from STDIN and writes statistics about the accuracy of the
* reconstruction to STDOUT
********************************************************************/

#include "../source/poisson_pedigree.h"
#include "../source/tree_diff_basic.h"
#include "../source/flags.h"

#include <iostream>

#define STOP_CHAR '~'

int main(int narg, char** args)
{

    // Construct pedigrees from STDIN
    std::string extant_dump, line;
    std::getline(std::cin, extant_dump, STOP_CHAR);
    poisson_pedigree* ped = poisson_pedigree::recover_dumped(extant_dump, new poisson_pedigree());
    std::getline(std::cin, extant_dump, STOP_CHAR);
    poisson_pedigree* rec = poisson_pedigree::recover_dumped(extant_dump, new poisson_pedigree());

    // Prepare the tree diff object
    tree_diff_basic* diff = new tree_diff_basic(ped, rec);
    diff->settings = LOG_WORK | LOG_DATA;

    // Flag definitions
    flag_reader fr;
    LOG_FLAG_READ(fr, diff);
    fr.add_flag("acc", 'a', 1, [&](std::vector<std::string> v, void* p) { diff->set_ch_acc(std::stod(v[0])); });
    if (fr.read_flags(narg, args) != FLAGS_INPUT_SUCCESS) {
        std::cout << "Invalid commands" << std::endl;
        return 1;
    }

    // Run tree diff and output
    diff->init();
    diff->topology_biject()->blocks_check();
    for (int i = 1; i < ped->num_grade(); i++)
        std::cout << "GENERATION " << i << ":\n" << tree_diff::stats_fmt(diff->nodes_total_gen[i], diff->nodes_correct_gen[i], diff->edges_total_gen[i], diff->edges_correct_gen[i], diff->blocks_total_gen[i], diff->blocks_attempted_gen[i], diff->blocks_correct_gen[i]) << std::endl;
    std::cout << "TOTAL:\n" << tree_diff::stats_fmt(diff->nodes_total, diff->nodes_correct, diff->edges_total, diff->edges_correct, diff->blocks_total, diff->blocks_attempted, diff->blocks_correct) << std::endl;
    return 0;

}
