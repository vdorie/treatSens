#ifndef CIBART_GUESS_NUM_CORES
#define CIBART_GUESS_NUM_CORES

#include <dbarts/cstdint.hpp>

namespace cibart {
  void guessNumCores(std::uint32_t* numPhyiscalProcessors, std::uint32_t* numLogicalProcessors);
}

#endif
