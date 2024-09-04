#include <vector>

#include "LifeAPI/LifeAPI.hpp"
#include "LifeAPI/Symmetry.hpp"

#include "Bloom.hpp"
#include "Countdown.hpp"
#include "Params.hpp"

const bool debug = false;
const bool print_progress = true;
const unsigned print_progress_frequency = 100000;

const unsigned approachRadius = 1; // Needs to match catalyst input
// const unsigned perturbationLookahead = 5; // TODO
const unsigned bloomPopulationThreshold = 12; // Min population

const unsigned maxStationaryGens = 32 - 1;

enum struct ContactType {
  CONTACT1,
  CONTACT2,
  CONTACTM,
  TRANSPARENT,
};

std::ostream &operator<<(std::ostream &out, const ContactType value) {
  return out << [value]() {
    switch (value) {
    case ContactType::CONTACT1:    return "CONTACT1";
    case ContactType::CONTACT2:    return "CONTACT2";
    case ContactType::CONTACTM:    return "CONTACTM";
    case ContactType::TRANSPARENT: return "TRANSPARENT";
    }
  }();
}

struct Approach {
  LifeState approachOn;
  LifeState approachOff;
  uint64_t signature;
  uint64_t signatureMask;

  void Move(int x, int y) {
    approachOn.Move(x, y);
    approachOff.Move(x, y);
    RecalculateSignature();
  }

  Approach Transformed(SymmetryTransform t) const {
    Approach result = {approachOn.Transformed(t),
                       approachOff.Transformed(t),
                       0,
                       0};
    result.RecalculateSignature();
    return result;
  }

  void RecalculateSignature() {
    signature = approachOn.GetPatch<approachRadius>({0, 0});
    signatureMask = (approachOn | approachOff).GetPatch<approachRadius>({0, 0});
  }

  bool MatchesSignature(uint64_t other) const {
    return ((signature ^ other) & signatureMask) == 0;
  }

  static Approach FromParams(const LifeHistoryState &params);
};

Approach Approach::FromParams(const LifeHistoryState &params) {
  LifeState catalyst = params.state & ~params.marked;

  LifeState reactionWithCatalyst = params.state;
  reactionWithCatalyst.Step();
  reactionWithCatalyst &= ~catalyst;

  LifeState reactionWithoutCatalyst = params.state & params.marked;
  reactionWithoutCatalyst.Step();

  Approach result;
  result.approachOn = params.marked & params.state;
  result.approachOff = (params.marked & ~params.state) | catalyst; // TODO:  | (result.required & ~result.state);

  assert((result.approachOn & result.approachOff).IsEmpty());

  return result;
}

struct CatalystData {
  LifeState state;
  LifeState halo;
  LifeState required;
  LifeState history1;
  LifeState history2;
  LifeState historyM;
  LifeState historyFlipped1;
  LifeState historyFlipped2;
  LifeState historyFlippedM;
  Approach approach;
  ContactType contactType; // At the origin contact cell
  unsigned maxRecoveryTime;
  unsigned minRecoveryTime;
  bool transparent;

  std::vector<Approach> forbidden;

  static CatalystData FromParamsNormal(CatalystParams &params);
  static CatalystData FromParamsTransparent(CatalystParams &params);
  static std::vector<CatalystData> FromParams(CatalystParams &params);

  CatalystData Transformed(SymmetryTransform t);

  LifeState CollisionMask(const CatalystData &b) const;
};

CatalystData CatalystData::FromParamsNormal(CatalystParams &params) {
  LifeHistoryState &approach = params.approaches[0];

  LifeState contact;
  {
    LifeState catalyst = approach.state & ~approach.marked;

    LifeState reactionWithCatalyst = approach.state;
    reactionWithCatalyst.Step();
    reactionWithCatalyst &= ~catalyst;

    LifeState reactionWithoutCatalyst =
        approach.state & approach.marked;
    reactionWithoutCatalyst.Step();

    contact = reactionWithCatalyst ^ reactionWithoutCatalyst;
  }

  auto contactorigin = contact.FirstOn();
  // TODO: Calculate approach from the soups?
  // This tramples the contents of `params`, seems bad
  approach.Move(-contactorigin.first, -contactorigin.second);

  params.state.AlignWith(approach.state & ~approach.marked);
  params.required.AlignWith(params.state);

  CatalystData result;

  result.state = params.state;
  result.halo = params.state.ZOI() & ~params.state;
  result.required = params.required.marked;
  result.state.InteractionCounts(result.history1, result.history2,
                                 result.historyM);
  result.historyFlipped1 = result.history1.Mirrored();
  result.historyFlipped2 = result.history2.Mirrored();
  result.historyFlippedM = result.historyM.Mirrored();
  result.approach.approachOn = approach.marked & approach.state;
  result.approach.approachOff = (approach.marked & ~approach.state) |
                       result.state | (result.required & ~result.state);

  assert((result.approach.approachOn & result.approach.approachOff).IsEmpty());

  for (auto &forbiddenpat : params.forbiddens) {
    Approach approach;
    forbiddenpat.AlignWith(params.state);
    approach.approachOn = forbiddenpat.state & forbiddenpat.marked;
    approach.approachOff = ~forbiddenpat.state & forbiddenpat.marked;
    approach.RecalculateSignature();
    result.forbidden.push_back(approach);
  }

  unsigned contactCount = result.approach.approachOn.CountNeighbours({0, 0});
  result.contactType =
    contactCount == 1
      ? ContactType::CONTACT1
      : (contactCount == 2 ? ContactType::CONTACT2 : ContactType::CONTACTM);

  result.approach.RecalculateSignature();

  result.minRecoveryTime = params.minRecoveryTime;
  result.maxRecoveryTime = params.maxRecoveryTime;

  result.transparent = false;

  if constexpr (debug) {
    std::cout << "Loaded catalyst: " << result.state << std::endl;
    // std::cout << "Params Approach: " << params.approach << std::endl;
    // std::cout << "Params Required: " << params.required << std::endl;
    std::cout << "Approach On: " << LifeHistoryState(result.state, LifeState(), result.approach.approachOn) << std::endl;
    std::cout << "Approach Off: " << LifeHistoryState(result.state, LifeState(), result.approach.approachOff) << std::endl;
    std::cout << "Required: " << LifeHistoryState(result.state, LifeState(), result.required) << std::endl;
    std::cout << "Contact Type: " << result.contactType << std::endl;
  }

  return result;
}

CatalystData CatalystData::FromParamsTransparent(CatalystParams &params) {
  CatalystData result;

  result.state = params.state;
  result.halo = params.state.ZOI() & ~params.state;
  result.required = LifeState();

  result.state.InteractionCounts(result.history1, result.history2,
                                 result.historyM);
  result.historyFlipped1 = result.history1.Mirrored();
  result.historyFlipped2 = result.history2.Mirrored();
  result.historyFlippedM = result.historyM.Mirrored();

  result.approach.approachOn = LifeState();
  result.approach.approachOff = LifeState();

  result.approach.signature = 0;
  result.approach.signatureMask = 0;

  result.contactType = ContactType::TRANSPARENT;

  result.minRecoveryTime = params.minRecoveryTime;
  result.maxRecoveryTime = params.maxRecoveryTime;

  result.transparent = true;

  if constexpr (debug) {
    std::cout << "Loaded transparent catalyst: " << result.state << std::endl;
  }

  return result;
}

CatalystData CatalystData::Transformed(SymmetryTransform t) {
  std::vector<Approach> newForbiddens;
  for (auto &oldForbidden : forbidden) {
    newForbiddens.push_back(oldForbidden.Transformed(t));
  }

  CatalystData result = {
    state.Transformed(t),
    halo.Transformed(t),
    required.Transformed(t),
    history1.Transformed(t),
    history2.Transformed(t),
    historyM.Transformed(t),
    historyFlipped1.Transformed(t),
    historyFlipped2.Transformed(t),
    historyFlippedM.Transformed(t),
    approach.Transformed(t),
    contactType,
    maxRecoveryTime,
    minRecoveryTime,
    transparent,
    newForbiddens
  };

  return result;
}

std::vector<CatalystData> CatalystData::FromParams(CatalystParams &params) {
  if (params.transparent) {
    CatalystData data = CatalystData::FromParamsTransparent(params);

    std::vector<CatalystData> result;
    for (auto t : params.state.SymmetryOrbitRepresentatives()) {
      result.push_back(data.Transformed(t));
    }
    return result;
  }

  CatalystData data = CatalystData::FromParamsNormal(params);
  std::vector<CatalystData> result;

  // TODO: use the actual symmetry of the catalyst
  using enum SymmetryTransform;
  for (auto t :
       {Identity, ReflectAcrossX, ReflectAcrossYeqX, ReflectAcrossY,
        ReflectAcrossYeqNegXP1, Rotate90, Rotate270, Rotate180OddBoth}) {
    result.push_back(data.Transformed(t));
  }

  return result;
}

LifeState CatalystData::CollisionMask(const CatalystData &b) const {
  // Block immediate births
  LifeState result = state.InteractionOffsets(b.state);

  // Block active interacting with required
  result |= (state.ZOI() & ~required).Convolve(b.required.Mirrored());

  // And vice versa
  result |= (b.state.ZOI() & ~b.required).Mirrored().Convolve(required);

  return result;
}

// Any precomputed data, constant at every node
struct SearchData {
  std::vector<CatalystData> catalysts;
  std::vector<LifeState> collisionMasks;
  LifeBloom *bloom;
};

struct Placement {
  std::pair<int, int> pos;
  unsigned catalystIx;
  unsigned gen;
};

// A set of catalyst placements
// mvrnote: name? Solution?
struct Configuration {
  LifeState state;
  LifeState catalysts;
  LifeState required;

  unsigned numCatalysts;
  unsigned numTransparent;

  unsigned lastInteraction;

  std::vector<Placement> placements;
  std::vector<LifeTarget> targets; // Pre-shifted catalysts

  Configuration()
      : state{}, catalysts{}, required{}, numCatalysts{0}, numTransparent{0},
        lastInteraction{0}, placements{}, targets{} {}
};

enum struct ProblemType {
  NONE,
  WINNER,
  REQUIRED,
  FILTER,
  UNRECOVERED,
  TOO_LONG,
  NO_REACTION,
  NOT_TRANSPARENT,
  STATIONARY,
  BLOOM_SEEN,
};

std::ostream &operator<<(std::ostream &out, const ProblemType value) {
  return out << [value]() {
    switch (value) {
    case ProblemType::NONE:            return "NONE";
    case ProblemType::WINNER:          return "WINNER";
    case ProblemType::REQUIRED:        return "REQUIRED";
    case ProblemType::FILTER:          return "FILTER";
    case ProblemType::UNRECOVERED:     return "UNRECOVERED";
    case ProblemType::TOO_LONG:        return "TOO_LONG";
    case ProblemType::NO_REACTION:     return "NO_REACTION";
    case ProblemType::NOT_TRANSPARENT: return "NOT_TRANSPARENT";
    case ProblemType::STATIONARY:      return "STATIONARY";
    case ProblemType::BLOOM_SEEN:      return "BLOOM_SEEN";
    }
  }();
}

struct Problem {
  std::pair<int, int> cell;
  unsigned gen;
  ProblemType type;

  LifeState LightCone(unsigned gen);
};

std::ostream &operator<<(std::ostream &out, const Problem value) {
  return out << value.type << " on gen " << value.gen << " at ("
             << value.cell.first << ", " << value.cell.second << ")";
}

// The state of a configuration after stepping
struct Lookahead {
  LifeState state;
  LifeState everActive;
  LifeCountdown<maxStationaryGens> stationaryCountdown;

  unsigned gen;
  bool hasInteracted;
  std::vector<unsigned> missingTime;
  std::vector<bool> catalystHasInteracted;
  unsigned startTime;
  unsigned recoveredTime;

  Lookahead()
      : state{}, everActive{}, gen{0}, hasInteracted{false}, missingTime{},
        catalystHasInteracted{}, recoveredTime{0} {}

  void Step(const Configuration &config);
  Problem Problem(const SearchParams &params, const SearchData &data,
                  const Configuration &config) const;

  std::pair<LifeState, bool> BloomKey(const Configuration &config) const;
};

void Lookahead::Step(const Configuration &config) {
  state.Step();
  gen++;
  everActive |= state ^ config.state;

  if (stationaryCountdown.n != 0) {
    LifeState stationary = state & ~config.catalysts;
    stationaryCountdown.TickOrReset(stationary);
  }

  // Why is there no easy way to iterate with index? My kingdom for a `zip`
  unsigned i = 0;
  bool allPresent = true;
  for (auto &t : config.targets) {
    if (state.Contains(t)) {
      missingTime[i] = 0;
    } else {
      missingTime[i]++;
      catalystHasInteracted[i] = true;
      if (!hasInteracted) {
        startTime = gen;
        hasInteracted = true;
      }
      allPresent = false;
    }
    i++;
  }

  recoveredTime = allPresent ? (recoveredTime + 1) : 0;
}

Problem Lookahead::Problem(const SearchParams &params, const SearchData &data,
                           const Configuration &config) const {
  {
    LifeState requiredViolations = config.required & (state ^ config.catalysts);
    std::pair<int, int> cell = requiredViolations.FirstOn();
    if (cell != std::make_pair(-1, -1))
      return {cell, gen, ProblemType::REQUIRED};
  }

  if (params.maxStationaryTime != 0) {
    LifeState stationaryViolations = stationaryCountdown.finished;
    std::pair<int, int> cell = stationaryViolations.FirstOn();
    if (cell != std::make_pair(-1, -1))
      return {cell, gen, ProblemType::STATIONARY};
  }

  {
    for (unsigned i = 0; i < config.numCatalysts; i++) {
      if (missingTime[i] >
          data.catalysts[config.placements[i].catalystIx].maxRecoveryTime) {
        const LifeTarget &target = config.targets[i];

        std::pair<int, int> cell = (target.wanted & ~state).FirstOn();
        if (cell.first == -1 && cell.second == -1)
          cell = (target.unwanted & state).FirstOn();

        return {cell, gen, ProblemType::UNRECOVERED};
      }
    }
  }

  {
    for (unsigned i = 0; i < config.numCatalysts; i++) {
      if (!data.catalysts[config.placements[i].catalystIx].transparent)
        continue;

      if (catalystHasInteracted[i] && missingTime[i] == 0) {
        LifeState nonTransparentCells = config.targets[i].wanted & ~everActive;

        std::pair<int, int> cell = nonTransparentCells.FirstOn();
        if (cell.first != -1 && cell.second != -1)
          return {cell, gen, ProblemType::NOT_TRANSPARENT};
      }
    }
  }

  {
    if (gen > params.maxFirstActiveGen && !hasInteracted)
      return {{-1, -1}, gen, ProblemType::NO_REACTION};
  }

  {
    if (gen > startTime + params.maxActiveWindowGens && hasInteracted)
      return {{-1, -1}, gen, ProblemType::TOO_LONG};
  }

  {
    if (gen > params.minFirstActiveGen && hasInteracted &&
        recoveredTime > params.minStableTime) {
      return {{-1, -1}, gen, ProblemType::WINNER};
    }
  }

  if (params.useBloomFilter) {
    auto [key, valid] = BloomKey(config);
    if(valid) {
      bool seen = data.bloom->Lookup(key);
      if (seen)
        return {{-1, -1}, gen, ProblemType::BLOOM_SEEN};
    }
  }

  return {{-1, -1}, gen, ProblemType::NONE};
}

std::pair<LifeState, bool> Lookahead::BloomKey(const Configuration &config) const {
  LifeState toHash = state & ~config.catalysts;
  bool valid = (gen > config.lastInteraction + 2) && (toHash.GetPop() > bloomPopulationThreshold);

  return {toHash, valid};
}

// Contact points that are close enough to the current problem to
// have an effect on it.
LifeState Problem::LightCone(unsigned currentgen) {
  switch (type) {
  case ProblemType::REQUIRED:
  case ProblemType::FILTER:
  case ProblemType::UNRECOVERED:
  case ProblemType::NOT_TRANSPARENT:
  case ProblemType::STATIONARY:
    return LifeState::NZOIAround(cell, gen - currentgen - 1);
  case ProblemType::WINNER:
  case ProblemType::NO_REACTION:
  case ProblemType::TOO_LONG:
    return ~LifeState();
  case ProblemType::NONE:
  case ProblemType::BLOOM_SEEN:
    __builtin_unreachable();
  }
}

// mvrnote: name?
// Arranged so 0 -> 1 is increasing information
struct CatalystConstraints {
  LifeState tried;
  LifeState knownPlaceable;
  LifeState knownUnplaceable;
};

struct SearchNode {
  Configuration config;
  Lookahead lookahead;

  LifeState history1;
  LifeState history2;
  LifeState historyM;

  std::vector<CatalystConstraints> constraints;

  SearchNode(const SearchParams &params, const SearchData &data);

  void BlockEarlyInteractions(const SearchParams &params,
                              const SearchData &data);

  void Step(const SearchParams &params, const SearchData &data);
};

Problem TryAdvance(const SearchParams &params, const SearchData &data,
                   const Configuration &config, SearchNode &search) {
  if constexpr (debug) std::cout << "Trying to advance: " << search.lookahead.state << std::endl;

  while (true) {
    Problem problem = search.lookahead.Problem(params, data, config);

    if (problem.type != ProblemType::NONE)
      return problem;

    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
        currentCountM(UNINITIALIZED);
    search.lookahead.state.InteractionCounts(currentCount1, currentCount2,
                                             currentCountM);
    LifeState newContactPoints = (currentCount1 & ~search.history1) |
                                 (currentCount2 & ~search.history2) |
                                 (currentCountM & ~search.historyM);

    bool needAdvance = true;

    if (!newContactPoints.IsEmpty()) {
      // See whether any catalysts have placements left at all
      for (unsigned i = 0; i < data.catalysts.size(); i++) {
        const CatalystData &catalyst = data.catalysts[i];
        LifeState remainingPlacements =
            newContactPoints & ~(search.constraints[i].knownUnplaceable |
                                 search.constraints[i].tried);

        switch (catalyst.contactType) {
        case ContactType::CONTACT1:
          remainingPlacements &= currentCount1;
          break;
        case ContactType::CONTACT2:
          remainingPlacements &= currentCount2;
          break;
        case ContactType::CONTACTM:
          remainingPlacements &= currentCountM;
          break;
        case ContactType::TRANSPARENT:
          break;
        }

        if (!remainingPlacements.IsEmpty()) {
          needAdvance = false;
          break;
        }
      }
    }

    if (!needAdvance)
      break;


    // TODO: reduce duplication
    if (params.useBloomFilter) {
      auto [key, valid] = search.lookahead.BloomKey(config);
      if(valid)
        data.bloom->Insert(key);
    }
    search.lookahead.Step(search.config);
    search.history1 |= currentCount1;
    search.history2 |= currentCount2;
    search.historyM |= currentCountM;

    if constexpr (debug) std::cout << "Advanced to " << search.lookahead.state << std::endl;
  }

  // We can't advance the state any more, but we still need to find the future
  // problem.
  Lookahead lookahead = search.lookahead;

  while (true) {
    lookahead.Step(config);

    if constexpr (debug) std::cout << "Lookahead to " << lookahead.state << std::endl;

    Problem problem = lookahead.Problem(params, data, config);

    if (problem.type != ProblemType::NONE)
      return problem;

    // TODO: reduce duplication
    if (params.useBloomFilter) {
      auto [key, valid] = lookahead.BloomKey(config);

      if(valid)
        data.bloom->Insert(key);
    }
  }
}

enum struct PlacementValidity {
  VALID,
  INVALID_CONTACT, // Wrong number of active neighbours
  FAILED_CONTACT,  // Invalid within `approachRadius` of the origin
  FAILED_ELSEWHERE // Invalid outside that, so will need to be re-checked
};

std::ostream &operator<<(std::ostream &out, const PlacementValidity value) {
  return out << [value]() {
    switch (value) {
    case PlacementValidity::VALID:            return "VALID";
    case PlacementValidity::INVALID_CONTACT:  return "INVALID_CONTACT";
    case PlacementValidity::FAILED_CONTACT:   return "FAILED_CONTACT";
    case PlacementValidity::FAILED_ELSEWHERE: return "FAILED_ELSEWHERE";
    }
  }();
}

PlacementValidity TestPlacement(const SearchData &data, SearchNode &search,
                                const LifeState &state, Placement p,
                                ContactType contactType, uint64_t signature,
                                const LifeState &historyCount1,
                                const LifeState &historyCount2,
                                const LifeState &currentCount1,
                                const LifeState &currentCount2) {
  const CatalystData &catalyst = data.catalysts[p.catalystIx];

  // Check that the type of catalyst contact matches the
  // active pattern neighbour count for an contact to occur
  if (catalyst.contactType != contactType)
    return PlacementValidity::INVALID_CONTACT;

  if (!catalyst.approach.MatchesSignature(signature))
    return PlacementValidity::FAILED_CONTACT;

  LifeState centered = state.Moved(-p.pos.first, -p.pos.second);

  LifeState mismatches = (catalyst.approach.approachOn & ~centered) |
                         (catalyst.approach.approachOff & centered);

  if (!mismatches.IsEmpty()) {
    constexpr LifeState originMask =
        LifeState::NZOIAround({0, 0}, approachRadius);
    if ((originMask & mismatches).IsEmpty()) {
      return PlacementValidity::FAILED_ELSEWHERE;
    } else {
      return PlacementValidity::FAILED_CONTACT;
    }
  }

  // TODO: should probably check M too...

  // Check whether this catalyst actually would have interacted in a previous
  // generation
  LifeState pastinteractions =
      (catalyst.history1 & historyCount2.Moved(-p.pos.first, -p.pos.second)) |
      (catalyst.history2 & historyCount1.Moved(-p.pos.first, -p.pos.second));
  if (!pastinteractions.IsEmpty()) {
    constexpr LifeState originMask =
        LifeState::NZOIAround({0, 0}, approachRadius);
    if ((originMask & pastinteractions).IsEmpty()) {
      return PlacementValidity::FAILED_ELSEWHERE;
    } else {
      return PlacementValidity::FAILED_CONTACT;
    }
  }

  // Check whether the catalyst's `required` is violated next generation
  LifeState immediatebirths = (catalyst.history1 & currentCount2.Moved(-p.pos.first, -p.pos.second)) |
                              (catalyst.history2 & currentCount1.Moved(-p.pos.first, -p.pos.second));
  immediatebirths &= (catalyst.required & ~catalyst.state);
  if (!immediatebirths.IsEmpty())
    return PlacementValidity::FAILED_ELSEWHERE;

  // Check the forbidden approaches
  for (auto &f : catalyst.forbidden) {
    LifeState differences = (f.approachOn & ~centered) |
                            (f.approachOff & centered);
    if (differences.IsEmpty())
      return PlacementValidity::FAILED_ELSEWHERE;
  }

  return PlacementValidity::VALID;
}

// This function is now too spaghetti
std::vector<Placement> CollectPlacements(const SearchParams &params,
                                         const SearchData &data,
                                         SearchNode &search, Problem &problem) {
  if constexpr (debug) std::cout << "Starting CollectPlacements" << std::endl;

  std::vector<Placement> result;

  LifeState current = search.lookahead.state;
  LifeState currentHistory1 = search.history1;
  LifeState currentHistory2 = search.history2;
  LifeState currentHistoryM = search.historyM;

  if constexpr (debug) std::cout << "Current: " << current << std::endl;
  if constexpr (debug) std::cout << "history1: " << currentHistory1 << std::endl;
  if constexpr (debug) std::cout << "history2: " << currentHistory2 << std::endl;
  if constexpr (debug) std::cout << "historyM: " << currentHistoryM << std::endl;

  LifeState somePlaceable = LifeState();
  for (unsigned i = 0; i < data.catalysts.size(); i++) {
    const CatalystData &catalyst = data.catalysts[i];
    somePlaceable |= ~(search.constraints[i].tried | search.constraints[i].knownUnplaceable);
  }

  bool advanceable = true;

  for (unsigned gen = search.lookahead.gen; gen < problem.gen; gen++) {
    if constexpr (debug) std::cout << "Gen " << gen << " state: " << current << std::endl;

    bool hasPlacement = false;

    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
        currentCountM(UNINITIALIZED);
    current.InteractionCounts(currentCount1, currentCount2, currentCountM);

    LifeState newContactPoints = (currentCount1 & ~currentHistory1) |
                                 (currentCount2 & ~currentHistory2) |
                                 (currentCountM & ~currentHistoryM);
    LifeState lightcone = problem.LightCone(gen);

    if(!advanceable)
      newContactPoints &= somePlaceable & lightcone;

    for (auto cell = newContactPoints.FirstOn(); cell != std::make_pair(-1, -1);
         newContactPoints.Erase(cell), cell = newContactPoints.FirstOn()) {

      bool inLightcone = lightcone.Get(cell);

      ContactType contactType =
          currentCount1.Get(cell)
              ? ContactType::CONTACT1
              : (currentCount2.Get(cell) ? ContactType::CONTACT2
                                         : ContactType::CONTACTM);

      uint64_t signature = current.GetPatch<approachRadius>(cell);

      // Handle ordinary catalysts
      for (unsigned i = 0; i < data.catalysts.size(); i++) {
        const CatalystData &catalyst = data.catalysts[i];

        // Fast checks
        if(!advanceable){
          if (catalyst.contactType != contactType)
            continue;

          if (!catalyst.approach.MatchesSignature(signature)) {
            search.constraints[i].knownUnplaceable.Set(cell);
            continue;
          }

          if (search.constraints[i].tried.Get(cell) ||
              search.constraints[i].knownUnplaceable.Get(cell)) {
            continue;
          }
        }

        Placement p = {cell, i, gen};

        if (catalyst.contactType == contactType &&
            search.constraints[i].knownPlaceable.Get(cell)) {
          if(inLightcone)
            result.push_back(p);
          hasPlacement = true;
          continue;
        }

        if (catalyst.contactType != contactType &&
            search.constraints[i].knownPlaceable.Get(cell)) {
          continue;
        }

        if (advanceable && search.constraints[i].knownUnplaceable.Get(cell)) {
          search.constraints[i].tried.Set(cell);
          continue;
        }

        if (search.constraints[i].tried.Get(cell))
          continue;

        PlacementValidity validity = TestPlacement(
            data, search, current, p, contactType, signature, currentHistory1,
            currentHistory2, currentCount1, currentCount2);

        switch (validity) {
        case PlacementValidity::VALID:
          search.constraints[i].knownPlaceable.Set(cell);
          if(inLightcone)
            result.push_back(p);
          hasPlacement = true;
          break;
        case PlacementValidity::FAILED_CONTACT:
          search.constraints[i].knownUnplaceable.Set(cell);
          if (advanceable)
            search.constraints[i].tried.Set(cell);
          break;
        case PlacementValidity::FAILED_ELSEWHERE:
        case PlacementValidity::INVALID_CONTACT:
          if (advanceable)
            search.constraints[i].tried.Set(cell);
          break;
        }
      }

      // Handle transparent catalysts
      if (search.config.numTransparent < params.maxTransparent) {
        LifeState newHistory1 = currentHistory1 & newContactPoints;
        LifeState newHistory2 = currentHistory2 & newContactPoints;
        for (unsigned i = 0; i < data.catalysts.size(); i++) {
          const CatalystData &catalyst = data.catalysts[i];
          if (catalyst.contactType != ContactType::TRANSPARENT)
            continue;

          LifeState newContactPoints =
              catalyst.historyFlipped1.Convolve(newHistory2) |
              catalyst.historyFlipped2.Convolve(newHistory1);
          newContactPoints &= ~search.constraints[i].tried;

          for (auto cell = newContactPoints.FirstOn();
               cell != std::make_pair(-1, -1); newContactPoints.Erase(cell),
                    cell = newContactPoints.FirstOn()) {
            Placement p = {cell, i, gen};
            if(inLightcone)
              result.push_back(p);
            hasPlacement = true;
          }
        }
      }
    }

    currentHistory1 |= currentCount1;
    currentHistory2 |= currentCount2;
    currentHistoryM |= currentCountM;

    if (advanceable && !hasPlacement) {
      // TODO: reduce duplication
      if (params.useBloomFilter) {
        auto [key, valid] = search.lookahead.BloomKey(search.config);
        if(valid)
          data.bloom->Insert(key);
      }

      search.lookahead.Step(search.config);
      current = search.lookahead.state;
      search.history1 |= currentCount1;
      search.history2 |= currentCount2;
      search.historyM |= currentCountM;

      if constexpr (debug) std::cout << "Advanced to " << search.lookahead.state << std::endl;
    } else {
      current.Step();
    }

    advanceable = advanceable && !hasPlacement;
  }

  return result;
}

void MakePlacement(const SearchParams &params, const SearchData &data,
                   SearchNode &search, const Placement &placement) {
  const CatalystData &catalyst = data.catalysts[placement.catalystIx];

  const LifeState catalystState = catalyst.state.Moved(placement.pos);

  search.config.numCatalysts++;
  if(catalyst.transparent) search.config.numTransparent++;
  search.config.lastInteraction = std::max(search.config.lastInteraction, placement.gen);
  search.config.state |= catalystState;
  search.config.catalysts |= catalystState;
  search.config.required |= catalyst.required.Moved(placement.pos);
  search.config.placements.push_back(placement);
  search.config.targets.push_back(
      LifeTarget(catalyst.state.Moved(placement.pos),
                 catalyst.halo.Moved(placement.pos)));

  search.lookahead.state |= catalystState;
  search.lookahead.missingTime.push_back(0);
  search.lookahead.catalystHasInteracted.push_back(false);

  search.history1 |= catalyst.history1.Moved(placement.pos);
  search.history2 |= catalyst.history2.Moved(placement.pos);
  search.historyM |= catalyst.historyM.Moved(placement.pos);

  for (unsigned t = 0; t < data.catalysts.size(); t++) {
    search.constraints[t].tried |=
        data.collisionMasks[placement.catalystIx * data.catalysts.size() + t]
            .Moved(placement.pos.first, placement.pos.second);
  }
}

void ResetLightcone(const SearchParams &params, const SearchData &data,
                    SearchNode &search, const Placement &placement) {

  if(debug) std::cout << "Resetting lightcone" << std::endl;

  LifeState current = search.lookahead.state;

  LifeState safeContacts1, safeContacts2, safeContactsM;

  // +1, because the first perturbation happens the generation after the placement
  for (unsigned gen = search.lookahead.gen; gen < placement.gen + 1; gen++) {
    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED), currentCountM(UNINITIALIZED);
    current.InteractionCounts(currentCount1, currentCount2, currentCountM);
    safeContacts1 |= currentCount1;
    safeContacts2 |= currentCount2;
    safeContactsM |= currentCountM;

    current.Step();
  }

  LifeState tooClose = LifeState::NZOIAround(placement.pos, approachRadius);

  // Until the universe is covered
  // TODO: this is a lot of generations, is there a sensible time to stop
  // sooner?
  for (int i = 2 * approachRadius + 1; i < 32; i += 2) {
    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED), currentCountM(UNINITIALIZED);
    current.InteractionCounts(currentCount1, currentCount2, currentCountM);
    safeContacts1 |= currentCount1 & ~tooClose;
    safeContacts2 |= currentCount2 & ~tooClose;
    safeContactsM |= currentCountM & ~tooClose;

    if(debug) std::cout << LifeHistoryState(current, LifeState(), tooClose) << std::endl;;

    current.Step();

    // Is it faster to just recompute the `NZOIAround`?
    tooClose = tooClose.ZOI();
  }

  // Now invalidate every placement that isn't safe
  // (Doesn't reset the `tried` field)
  for (unsigned i = 0; i < data.catalysts.size(); i++) {
    switch (data.catalysts[i].contactType) {
    case ContactType::CONTACT1:
      search.constraints[i].knownPlaceable &= safeContacts1;
      search.constraints[i].knownUnplaceable &= safeContacts1;
      break;
    case ContactType::CONTACT2:
      search.constraints[i].knownPlaceable &= safeContacts2;
      search.constraints[i].knownUnplaceable &= safeContacts2;
      break;
    case ContactType::CONTACTM:
      search.constraints[i].knownPlaceable &= safeContactsM;
      search.constraints[i].knownUnplaceable &= safeContactsM;
      break;
    case ContactType::TRANSPARENT:
      // They are always placeable
      break;
    }
  }
}

void RunSearch(const SearchParams &params, const SearchData &data,
               SearchNode &search) {

  if constexpr (debug) std::cout << "Starting node: " << search.config.state << std::endl;

  Problem problem = TryAdvance(params, data, search.config, search);

  if constexpr (print_progress) {
    static unsigned counter = 0;
    counter++;
    if (counter == print_progress_frequency) [[unlikely]] {
      std::cout << "Current configuration: " << search.config.state << std::endl;
      LifeState problemGen = search.lookahead.state;
      problemGen.Step(problem.gen - search.lookahead.gen);
      std::cout << "Current problem: " << problem << std::endl;
      if(problem.cell.first != -1)
        std::cout << LifeHistoryState(problemGen, LifeState(), LifeState::Cell(problem.cell)) << std::endl;
      if(params.useBloomFilter) {
        std::cout << "Bloom filter population: " << data.bloom->items << std::endl;
        std::cout << "Bloom filter error rate: " << data.bloom->ApproximateErrorRate() << std::endl;
      }
      std::cout << "Current placements:" << std::endl;
      LifeState progression = params.state.state;
      for (auto &p : search.config.placements) {
        const CatalystData &catalystdata = data.catalysts[p.catalystIx];
        const LifeState catalyst = catalystdata.state.Moved(p.pos);
        progression |= catalyst;
        std::cout << progression << std::endl;
      }
      counter = 0;
    }
  }

  if constexpr (debug) {
    std::cout << "Problem: " << problem << std::endl;
    if (problem.cell.first != -1) {
      LifeState problemGen = search.lookahead.state;
      problemGen.Step(problem.gen - search.lookahead.gen);
      std::cout << LifeHistoryState(problemGen, LifeState(), LifeState::Cell(problem.cell)) << std::endl;
    }
  }

  if (problem.type == ProblemType::BLOOM_SEEN)
    return;
  
  if (problem.type == ProblemType::TOO_LONG) {
    std::cout << "Too long: " << search.config.state << std::endl;
  }

  if (problem.type == ProblemType::WINNER) {
    // TODO report properly!
    std::cout << "Winner: " << search.config.state << std::endl;
    if constexpr (debug) {
      std::cout << "Required: " << LifeHistoryState(search.config.state, LifeState(), search.config.required) << std::endl;
      for (auto &p : search.config.placements) {
        std::cout << "Placement: " << p.catalystIx << " at (" << p.pos.first << ", " << p.pos.second << ")" << std::endl;
      }
    }
  }

  if (search.config.numCatalysts == params.maxCatalysts)
    return;

  std::vector<Placement> placements =
      CollectPlacements(params, data, search, problem);

  for (const auto &placement : placements) {
    // TODO: there may be repeat placements due to transparent catalysts
    if (search.constraints[placement.catalystIx].tried.Get(placement.pos))
      continue;

    search.constraints[placement.catalystIx].tried.Set(placement.pos);

    SearchNode newSearch = search;

    MakePlacement(params, data, newSearch, placement);

    if constexpr (debug) {
      if (params.hasOracle && !(newSearch.config.state & ~params.oracle).IsEmpty()) {
        if constexpr (debug)
          std::cout << "Oracle failed: " << newSearch.lookahead.state << std::endl;
        continue;
      }
    }

    if constexpr (debug)
      std::cout << "Branching: " << placement.catalystIx << std::endl;

    ResetLightcone(params, data, newSearch, placement);

    RunSearch(params, data, newSearch);
  }
}

std::vector<LifeState>
CalculateCollisionMasks(const std::vector<CatalystData> &catalysts) {
  unsigned count = catalysts.size();
  std::vector<LifeState> result(count * count);
  for (unsigned s = 0; s < count; s++) {
    for (unsigned t = 0; t < count; t++) {
      result[s * count + t] = catalysts[s].CollisionMask(catalysts[t]);
    }
  }
  return result;
}

void SearchNode::BlockEarlyInteractions(const SearchParams &params,
                                        const SearchData &data) {
  LifeState statezoi = params.state.state.ZOI();

  // Always block overlapping placements
  for (unsigned i = 0; i < data.catalysts.size(); i++) {
    const CatalystData &catalyst = data.catalysts[i];
    constraints[i].tried |= catalyst.state.Mirrored().Convolve(statezoi);
  }

  LifeState current = params.state.state;
  // Also block early interactions
  for (unsigned g = 0; g < params.minFirstActiveGen; g++) {
    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
        currentCountM(UNINITIALIZED);
    current.InteractionCounts(currentCount1, currentCount2, currentCountM);
    for (unsigned i = 0; i < data.catalysts.size(); i++) {
      const CatalystData &catalyst = data.catalysts[i];
      constraints[i].tried |= catalyst.historyFlipped1.Convolve(currentCount2) |
                              catalyst.historyFlipped2.Convolve(currentCount1);
    }
    current.Step();
  }
}

SearchNode::SearchNode(const SearchParams &params, const SearchData &data) {
  config = Configuration();
  lookahead = Lookahead();

  config.state = params.state.state;
  lookahead.state = params.state.state;
  lookahead.stationaryCountdown.n = params.maxStationaryTime;

  config.state.InteractionCounts(history1, history2, historyM);
  constraints = std::vector<CatalystConstraints>(data.catalysts.size(),
                                                 CatalystConstraints());

  if (!params.state.history.IsEmpty()) {
    for (unsigned i = 0; i < data.catalysts.size(); i++) {
      constraints[i].tried |=
          data.catalysts[i].state.Mirrored().Convolve(~params.state.history);
    }
  }
}

int main(int, char *argv[]) {
  auto toml = toml::parse(argv[1]);
  SearchParams params = SearchParams::FromToml(toml);

  std::vector<CatalystData> catalystdata;
  for (auto &c : params.catalysts) {
    auto newdata = CatalystData::FromParams(c);
    catalystdata.insert(catalystdata.end(), newdata.begin(), newdata.end()); // Why is C++ like this
  }

  std::vector<LifeState> masks = CalculateCollisionMasks(catalystdata);

  LifeBloom *bloom = 0;
  if (params.useBloomFilter)
    bloom = new LifeBloom();

  SearchData data = {catalystdata, masks, bloom};

  SearchNode search(params, data);

  search.BlockEarlyInteractions(params, data);

  RunSearch(params, data, search);
  if (params.useBloomFilter) {
    std::cout << "Bloom filter population: " << data.bloom->items << std::endl;
    std::cout << "Bloom filter approx    : "
              << data.bloom->ApproximatePopulation() << std::endl;
    std::cout << "Bloom filter error rate: "
              << data.bloom->ApproximateErrorRate() << std::endl;
  }
}
