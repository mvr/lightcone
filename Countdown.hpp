#pragma once

#include "LifeAPI/LifeAPI.hpp"

// Copied from Barrister
template <uint32_t max>
class LifeCountdown {
public:
  static constexpr uint32_t lmax = (max == 0) ? 0 : 32 - __builtin_clz(max);
  LifeState started;
  LifeState finished;
  std::array<LifeState, lmax> counter = {0};
  uint32_t n;

  LifeCountdown() : started{}, finished{}, counter{}, n{0} {};
  LifeCountdown(uint32_t n) : started{}, finished{}, counter{}, n{n} {};

  void Start(const LifeState &state) {
    LifeState newStarted = state & ~started;
    for (unsigned i = 0; i < lmax; i++) {
      if ((n >> i) & 1) {
        counter[i] |= newStarted;
      }
    }
    started |= state;
  }

  void Reset(const LifeState &state) {
    for (unsigned i = 0; i < lmax; i++) {
      counter[i] &= ~state;
    }
    started &= ~state;
  }

  void Tick() {
    auto carry = started;
    for (unsigned i = 0; i < lmax; i++) {
        counter[i] ^= carry;
        carry &= counter[i];
    }
    finished |= carry;
  }

  void TickOrReset(const LifeState &state) {
    LifeState newStarted = state & ~started;
    LifeState carry = state & started;

    for (unsigned i = 0; i < lmax; i++) {
      if ((n >> i) & 1) {
        counter[i] |= newStarted;
      } else {
        counter[i] &= ~newStarted;
      }
      counter[i] ^= carry;
      carry &= counter[i];
    }
    finished |= carry;
    started = state;
  }
};
