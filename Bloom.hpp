#pragma once

#include <bit>
#include <cmath>

#include "assert.h"

#include "LifeAPI/LifeAPI.hpp"

const unsigned bloomHashes = 12;
const unsigned bloomSize = 8 * 1024 * 1024 * 16; // 16MiB
const uint64_t bloomSeed = 19;

struct LifeBloom {
  uint64_t table[bloomSize / 64];

  unsigned items;

  LifeBloom() : table{0}, items{0} {}
  void Insert(const LifeState &state);
  bool Lookup(const LifeState &state) const;
  bool InsertAndLookup(const LifeState &state);

  unsigned ApproximatePopulation() const;
  float ApproximateErrorRate() const;
};

inline void LifeBloom::Insert(const LifeState &state) {
  auto [baseHash1, baseHash2] = XXH3_128bits_withSeed(state.state, 8 * N, bloomSeed);

  for (unsigned i = 0; i < bloomHashes; i++) {
    unsigned hash = (baseHash1 + i * baseHash2) % bloomSize;
    table[hash / 64] |= 1ULL << (hash % 64);
  }

  items++;
}

inline bool LifeBloom::Lookup(const LifeState &state) const {
  auto [baseHash1, baseHash2] = XXH3_128bits_withSeed(state.state, 8 * N, bloomSeed);

  for (unsigned i = 0; i < bloomHashes; i++) {
    unsigned hash = (baseHash1 + i * baseHash2) % bloomSize;
    bool bit = table[hash / 64] & (1ULL << (hash % 64));
    if (bit == 0)
      return false;
  }
  return true;
}

inline bool LifeBloom::InsertAndLookup(const LifeState &state) {
  auto [baseHash1, baseHash2] = XXH3_128bits_withSeed(state.state, 8 * N, bloomSeed);

  bool seen = true;
  for (unsigned i = 0; i < bloomHashes; i++) {
    unsigned hash = (baseHash1 + i * baseHash2) % bloomSize;
    uint64_t bit = table[hash / 64] & (1ULL << (hash % 64));
    if (bit == 0)
      seen = false;
    table[hash / 64] |= 1ULL << (hash % 64);
  }
  return seen;
}

inline unsigned LifeBloom::ApproximatePopulation() const {
  unsigned setBits = 0;
  for (unsigned i = 0; i < bloomSize / 64; i++) {
    setBits += std::popcount(table[i]);
  }
  float approx = -(float)bloomSize / (float)bloomHashes *
                 std::log1p(-(float)setBits / (float)bloomSize);

  return (unsigned)approx;
}

inline float LifeBloom::ApproximateErrorRate() const {
  float perhash = -std::expm1(-(float)bloomHashes * (float)items / (float)bloomSize);
  return std::pow(perhash, bloomHashes);
}
