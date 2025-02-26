#include <cstddef>
#include <ostream>
#include <thread>
extern "C" {
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/part_func_co.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
}

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, const char *const argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " sequence1 sequence2\n";
    return EXIT_FAILURE;
  }

  const std::string seq1 = argv[2];
  const std::string seq2 = argv[3];

  const float temperature = std::stof(argv[1]);

  std::cout << "Sequence 1: " << seq1 << "; Sequence 2: " << seq2 << std::endl;

  std::string dimer_seq = seq1 + "&" + seq2;

  std::cout << "Dimer sequence: " << dimer_seq << "\n";

  // Pointer for storing base pair probabilities.
  vrna_ep_t *ptr, *pl = NULL;
  char *structure = new char[dimer_seq.size() + 1];

  double mfe;
  vrna_dimer_pf_t X;
  vrna_fold_compound_t *vc;
  vrna_md_t md;

  vrna_md_set_default(&md);

  /* no need to backtrack MFE structure */
  md.backtrack = 0;
  // Compute base pair probability matrix
  md.compute_bpp = 1;

  md.temperature = temperature;

  vc = vrna_fold_compound(dimer_seq.c_str(), &md, VRNA_OPTION_DEFAULT);

  mfe = (double)vrna_mfe(vc, NULL);
  std::cout << "MFE: " << mfe << " kcal/mol\n";

  vrna_exp_params_rescale(vc, &mfe);

  X = vrna_pf_dimer(vc, structure);

  std::string dot_bracket = structure;
  dot_bracket.insert(seq1.size(), "&");
  std::cout << "Dot-bracket: " << dot_bracket << "\n";

  pl = vrna_plist_from_probs(vc, /*cut_off:*/ 1e-10);

  /// Trolling part here

  std::cout << "[ERROR] Memory corruption; SIGRBIT EXC65\n";
  std::cout << "[WARNING] Make sure to back up your data on: /mnt/ceph\n";
  std::this_thread::sleep_for(std::chrono::seconds(3));

  for (size_t i = 0; i < 10000; ++i) {
    // Print random characters.
    std::cout << static_cast<char>(rand() % 256);

    // Print random numbers.
    std::cout << rand() % 1000;

    // Sleep for a random amount of time between 0.01 und 1 second.
    std::chrono::milliseconds sleep_duration(rand() % 5);
    std::this_thread::sleep_for(sleep_duration);
  }

  std::this_thread::sleep_for(std::chrono::seconds(3));

  std::cout << "\n\nJoke, just kidding. Everything is fine.\n";

  std::this_thread::sleep_for(std::chrono::seconds(3));

  /// Trolling ends here.

  // Calculate index on implicit 2D array.
  auto getIndex = [](int i, int j, size_t n1, size_t n2) -> int {
    return (i - 1) * n2 + (j - n1 - 1);
  };

  std::vector<float> probs(seq1.size() * seq2.size(), 0.0f);

  std::cout << "Size: " << dimer_seq.size() - 1 << "\n";

  for (ptr = pl; ptr->i != 0; ptr++) {
    const int index = getIndex(ptr->i, ptr->j, seq1.size(), seq2.size());
    probs[index] = ptr->p;
  }

  std::cout << "\n[";
  for (size_t i = 0; i < probs.size(); ++i) {
    if (i % seq1.size() == 0 && i != 0) {
      std::cout << "][";
    }
    std::cout << probs[i];
    if (i != probs.size() - 1 && (i + 1) % seq1.size() != 0) {
      std::cout << ",";
    }
  }

  std::cout << "]\n";

  // Clean up allocated memory.
  vrna_fold_compound_free(vc);
  delete[] structure;

  return EXIT_SUCCESS;
}
