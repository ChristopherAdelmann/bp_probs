#include <cstring>
extern "C" {
#include <ViennaRNA/cofold.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/subopt.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
}
#include <cstdlib>
#include <iostream>

auto main(int argc, const char *const argv[]) -> int {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " sequence1&sequence2\n";
        return EXIT_FAILURE;
    }

    const char *seq = argv[1];
    std::cout << "Sequence: " << seq << std::endl;

    // Allocate a buffer for the structure string.
    size_t seq_length = std::strlen(seq);
    char *structure = new char[seq_length + 1];
    if (!structure) {
        std::cerr << "Allocation failed for structure buffer\n";
        return EXIT_FAILURE;
    }

    // Pointer for storing base pair probabilities.
    vrna_ep_t *pair_probabilities = nullptr;

    // Compute the partition function for the co-folded sequence.
    auto energy = vrna_pf_co_fold(seq, structure, &pair_probabilities);

    // std::cout << structure << " [ " << energy << " ]" << std::endl;

    // Optionally, print base pairs with high probability.
    for (vrna_ep_t *ptr = pair_probabilities; ptr && ptr->i != 0; ++ptr) {
        std::cout << "(" << ptr->i << ", " << ptr->j << ") : " << ptr->p << "\n";
    }

    // Clean up allocated memory.
    delete[] structure;
    free(pair_probabilities);

    return EXIT_SUCCESS;
}
