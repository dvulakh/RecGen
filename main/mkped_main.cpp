
/********************************************************************
* Stochastically generates a poisson pedigree based on properties
* read from command-line arguments and prints it to STDOUT.
********************************************************************/

#include "../source/poisson_pedigree.h"
#include "../source/flags.h"

#include <iostream>
#include <iterator>

#define STOP_CHAR '~'

int main(int narg, char** args)
{

    // If there are no arguments, parse the shorthand
    if (narg == 1) {
        std::cin >> std::noskipws;
        std::istream_iterator<char> it(std::cin), end;
        poisson_pedigree* ped = poisson_pedigree::parse_shorthand(std::string(it, end));
        std::cout << ped->dump_extant() << std::endl << STOP_CHAR << std::endl << ped->dump() << std::endl;
        return 0;
    }

    // Read pedigree properties
    std::string arg;
    for (int i = 1; i < narg; i++)
        arg += std::string(args[i]) + " ";
    poisson_pedigree* ped = poisson_pedigree::recover_dumped(arg, new poisson_pedigree());

    // Generate and print pedigree
    std::cout << ped->build()->dump_extant() << std::endl << STOP_CHAR << std::endl << ped->dump() << std::endl;
    return 0;

}
