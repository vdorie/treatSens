#ifndef CIBART_GUESS_NUM_CORES
#define CIBART_GUESS_NUM_CORES

#include <cstdint>

namespace cibart {
  void guessNumCores(std::uint32_t* numPhyiscalProcessors, std::uint32_t* numLogicalProcessors);
}

#endif
