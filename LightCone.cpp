#include <vector>

#include "LifeAPI/LifeAPI.hpp"
#include "LifeAPI/Symmetry.hpp"

#include "Bloom.hpp"
#include "Countdown.hpp"
#include "Params.hpp"

const bool debug = false;
const bool print_progress = true;

const unsigned approachRadius = 1; // Needs to match catalyst input
// const unsigned perturbationLookahead = 5; // TODO
const unsigned bloomThreshold = 12; // Min population

const unsigned maxStationaryGens = 32 - 1;

enum ContactType {
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

struct CatalystData {
  LifeState state;
  LifeState halo;
  LifeState required;
  LifeState history1;
  LifeState history2;
  LifeState historyM;
  LifeState approachOn;
  LifeState approachOff;
  LifeState historyFlipped1;
  LifeState historyFlipped2;
  LifeState historyFlippedM;
  ContactType contactType; // At the origin contact cell
  uint64_t signature;
  uint64_t signatureMask;
  unsigned maxRecoveryTime;
  unsigned minRecoveryTime;
  bool transparent;

  static CatalystData FromParamsNormal(CatalystParams &params);
  static CatalystData FromParamsTransparent(CatalystParams &params);
  static std::vector<CatalystData> FromParams(CatalystParams &params);

  CatalystData Transformed(SymmetryTransform t);

  LifeState CollisionMask(const CatalystData &b) const;

  bool MatchesSignature(uint64_t other) const {
    return ((signature ^ other) & signatureMask) == 0;
  }
};

CatalystData CatalystData::FromParamsNormal(CatalystParams &params) {
  LifeState contact;
  {
    LifeState catalyst = params.approach.state & ~params.approach.marked;

    LifeState reactionWithCatalyst = params.approach.state;
    reactionWithCatalyst.Step();
    reactionWithCatalyst &= ~catalyst;

    LifeState reactionWithoutCatalyst =
        params.approach.state & params.approach.marked;
    reactionWithoutCatalyst.Step();

    contact = reactionWithCatalyst ^ reactionWithoutCatalyst;
  }

  auto contactorigin = contact.FirstOn();
  // TODO: Calculate approach from the soups?
  // This tramples the contents of `params`, seems bad
  params.approach.Move(-contactorigin.first, -contactorigin.second);

  params.state.AlignWith(params.approach.state & ~params.approach.marked);
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
  result.approachOn = params.approach.marked & params.approach.state;
  result.approachOff = (params.approach.marked & ~params.approach.state) |
                       result.state | (result.required & ~result.state);

  assert((result.approachOn & result.approachOff).IsEmpty());

  unsigned contactCount = result.approachOn.CountNeighbours({0, 0});
  result.contactType =
      contactCount == 1
          ? ContactType::CONTACT1
          : (contactCount == 2 ? ContactType::CONTACT2 : ContactType::CONTACTM);

  result.signature = result.approachOn.GetPatch<approachRadius>({0, 0});
  result.signatureMask =
      (result.approachOn | result.approachOff).GetPatch<approachRadius>({0, 0});

  result.minRecoveryTime = params.minRecoveryTime;
  result.maxRecoveryTime = params.maxRecoveryTime;

  result.transparent = false;

  if constexpr (debug) {
    std::cout << "Loaded catalyst: " << result.state << std::endl;
    // std::cout << "Params Approach: " << params.approach << std::endl;
    // std::cout << "Params Required: " << params.required << std::endl;
    std::cout << "Approach On: " << LifeHistoryState(result.state, LifeState(), result.approachOn) << std::endl;
    std::cout << "Approach Off: " << LifeHistoryState(result.state, LifeState(), result.approachOff) << std::endl;
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

  result.approachOn = LifeState();
  result.approachOff = LifeState();

  result.contactType = ContactType::TRANSPARENT;

  result.signature = 0;
  result.signatureMask = 0;

  result.minRecoveryTime = params.minRecoveryTime;
  result.maxRecoveryTime = params.maxRecoveryTime;

  result.transparent = true;

  if constexpr (debug) {
    std::cout << "Loaded transparent catalyst: " << result.state << std::endl;
  }

  return result;
}

CatalystData CatalystData::Transformed(SymmetryTransform t) {
  CatalystData result = {
    state.Transformed(t),
    halo.Transformed(t),
    required.Transformed(t),
    history1.Transformed(t),
    history2.Transformed(t),
    historyM.Transformed(t),
    approachOn.Transformed(t),
    approachOff.Transformed(t),
    historyFlipped1.Transformed(t),
    historyFlipped2.Transformed(t),
    historyFlippedM.Transformed(t),
    contactType,
    signature,
    signatureMask,
    maxRecoveryTime,
    minRecoveryTime,
    transparent
  };

  // Needs to be recalculated
  result.signature = result.approachOn.GetPatch<approachRadius>({0, 0});
  result.signatureMask =
      (result.approachOn | result.approachOff).GetPatch<approachRadius>({0, 0});

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
  return state.InteractionOffsets(b.state);
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

  bool Bloomable(const Configuration &config) const;

  LifeState BloomKey(const Configuration &config) const;
};

void Lookahead::Step(const Configuration &config) {
  state.Step();
  gen++;
  everActive |= state ^ config.state;

  if (stationaryCountdown.n != 0) {
    LifeState stationary = state & ~config.catalysts;
    stationaryCountdown.Reset(~stationary);
    stationaryCountdown.Start(stationary);
    stationaryCountdown.Tick();
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

  if (params.useBloomFilter && Bloomable(config)) {
    LifeState key = BloomKey(config);
    if (key.GetPop() > bloomThreshold) {
      bool seen = data.bloom->Lookup(key);
      if (seen)
        return {{-1, -1}, gen, ProblemType::BLOOM_SEEN};
    }
  }

  return {{-1, -1}, gen, ProblemType::NONE};
}

bool Lookahead::Bloomable(const Configuration &config) const {
  return gen > config.lastInteraction + 2;
}

// TODO: this doesn't work well for catalysts with long internal recoveries
LifeState Lookahead::BloomKey(const Configuration &config) const {
  return state & ~config.catalysts;
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
    return LifeState::NZOIAround(cell, gen - currentgen);
  case ProblemType::WINNER:
  case ProblemType::NO_REACTION:
  case ProblemType::TOO_LONG:
  case ProblemType::BLOOM_SEEN:
    return ~LifeState();
  case ProblemType::NONE:
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

// TODO: determine exactly when this needs to be run.
// Is it after the final perturbation?
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
        if (!remainingPlacements.IsEmpty()) {
          needAdvance = false;
          break;
        }
      }
    }

    if (needAdvance) {
      search.lookahead.Step(search.config);
      search.history1 |= currentCount1;
      search.history2 |= currentCount2;
      search.historyM |= currentCountM;

      if constexpr (debug) std::cout << "Advanced to " << search.lookahead.state << std::endl;

      // TODO: reduce duplication
      if (params.useBloomFilter && search.lookahead.Bloomable(search.config)) {
        LifeState key = search.lookahead.BloomKey(search.config);
        if (key.GetPop() > bloomThreshold)
          data.bloom->Insert(key);
      }
    } else {
      break;
    }
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
    if (params.useBloomFilter && lookahead.Bloomable(search.config)) {
      LifeState key = lookahead.BloomKey(search.config);
      if (key.GetPop() > bloomThreshold)
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

  if (!catalyst.MatchesSignature(signature))
    return PlacementValidity::FAILED_CONTACT;

  LifeState mismatches = (catalyst.approachOn.Moved(p.pos) & ~state) |
                         (catalyst.approachOff.Moved(p.pos) & state);

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
      (catalyst.history1.Moved(p.pos) & historyCount2) |
      (catalyst.history2.Moved(p.pos) & historyCount1);
  if (!pastinteractions.IsEmpty()) {
    constexpr LifeState originMask =
        LifeState::NZOIAround({0, 0}, approachRadius);
    if ((originMask & pastinteractions).IsEmpty()) {
      return PlacementValidity::FAILED_ELSEWHERE;
    } else {
      return PlacementValidity::FAILED_CONTACT;
    }
  }

  // Check whether the cataylst's `required` is violated next generation
  LifeState immediatebirths = (catalyst.history1.Moved(p.pos) & currentCount2) |
                              (catalyst.history2.Moved(p.pos) & currentCount1);
  immediatebirths &= (catalyst.required & ~catalyst.state).Moved(p.pos);
  if (!immediatebirths.IsEmpty())
    return PlacementValidity::FAILED_ELSEWHERE;

  return PlacementValidity::VALID;
}

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

  for (unsigned g = search.lookahead.gen; g < problem.gen; g++) {
    if constexpr (debug) std::cout << "Gen " << g << " state: " << current << std::endl;

    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
        currentCountM(UNINITIALIZED);
    current.InteractionCounts(currentCount1, currentCount2, currentCountM);

    LifeState newContactPoints = (currentCount1 & ~currentHistory1) |
                                 (currentCount2 & ~currentHistory2) |
                                 (currentCountM & ~currentHistoryM);

    newContactPoints &= problem.LightCone(g) & somePlaceable;

    for (auto cell = newContactPoints.FirstOn(); cell != std::make_pair(-1, -1);
         newContactPoints.Erase(cell), cell = newContactPoints.FirstOn()) {

      ContactType contactType =
          currentCount1.Get(cell)
              ? CONTACT1
              : (currentCount2.Get(cell) ? CONTACT2 : CONTACTM);

      uint64_t signature = current.GetPatch<approachRadius>(cell);

      // Handle ordinary catalysts
      for (unsigned i = 0; i < data.catalysts.size(); i++) {
        const CatalystData &catalyst = data.catalysts[i];

        if (catalyst.contactType != contactType)
          continue;

        if (!catalyst.MatchesSignature(signature)) {
          search.constraints[i].knownUnplaceable.Set(cell);
          continue;
        }

        if (search.constraints[i].tried.Get(cell) ||
            search.constraints[i].knownUnplaceable.Get(cell)) {
          continue;
        }

        Placement p = {cell, i, g};

        if (catalyst.contactType == contactType &&
            search.constraints[i].knownPlaceable.Get(cell)) {
          result.push_back(p);
          continue;
        }

        PlacementValidity validity = TestPlacement(
            data, search, current, p, contactType, signature, currentHistory1,
            currentHistory2, currentCount1, currentCount2);

        switch (validity) {
        case PlacementValidity::VALID:
          search.constraints[i].knownPlaceable.Set(cell);
          result.push_back(p);
          break;
        case PlacementValidity::FAILED_CONTACT:
          search.constraints[i].knownUnplaceable.Set(cell);
          break;
        case PlacementValidity::FAILED_ELSEWHERE:
        case PlacementValidity::INVALID_CONTACT:
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
            Placement p = {cell, i, g};
            result.push_back(p);
          }
        }
      }
    }

    currentHistory1 |= currentCount1;
    currentHistory2 |= currentCount2;
    currentHistoryM |= currentCountM;

    current.Step();
  }

  return result;
}

void MakePlacement(const SearchParams &params, const SearchData &data,
                   SearchNode &search, const Placement &placement) {
  const CatalystData &catalystdata = data.catalysts[placement.catalystIx];

  const LifeState catalyst = catalystdata.state.Moved(placement.pos);

  search.config.numCatalysts++;
  if(catalystdata.transparent) search.config.numTransparent++;
  search.config.lastInteraction = std::max(search.config.lastInteraction, placement.gen);
  search.config.state |= catalyst;
  search.config.catalysts |= catalyst;
  search.config.required |= catalystdata.required.Moved(placement.pos);
  search.config.placements.push_back(placement);
  search.config.targets.push_back(
      LifeTarget(catalystdata.state.Moved(placement.pos),
                 catalystdata.halo.Moved(placement.pos)));

  search.lookahead.state |= catalyst;
  search.lookahead.missingTime.push_back(0);
  search.lookahead.catalystHasInteracted.push_back(false);

  search.history1 |= catalystdata.history1.Moved(placement.pos);
  search.history2 |= catalystdata.history2.Moved(placement.pos);
  search.historyM |= catalystdata.historyM.Moved(placement.pos);

  for (unsigned t = 0; t < data.catalysts.size(); t++) {
    search.constraints[t].tried |=
        data.collisionMasks[placement.catalystIx * data.catalysts.size() + t]
            .Moved(placement.pos.first, placement.pos.second);
  }
}

void ResetLightcone(const SearchParams &params, const SearchData &data,
                    SearchNode &search, const Placement &placement) {
  // TODO: Assuming that the placement has already been made
  LifeState current = search.lookahead.state;

  LifeState tooClose = LifeState::NZOIAround(placement.pos, approachRadius);

  LifeState safeContacts;

  // Until the universe is covered
  // TODO: this is a lot of generations, is there a sensible time to stop
  // sooner?
  for (int i = 2 * approachRadius + 1; i < 32; i += 2) {
    safeContacts |= current.ZOI() & ~tooClose;
    current.Step();
    // Is it faster to just recompute the `NZOIAround`?
    tooClose = tooClose.ZOI();
  }

  // Now invalidate every placement that isn't safe
  // (Doesn't reset the `tried` field)
  for (unsigned i = 0; i < data.catalysts.size(); i++) {
    search.constraints[i].knownPlaceable &= safeContacts;
    search.constraints[i].knownUnplaceable &= safeContacts;
  }
}

void RunSearch(const SearchParams &params, const SearchData &data,
               SearchNode &search) {

  if constexpr (debug) std::cout << "Starting node: " << search.lookahead.state << std::endl;

  Problem problem = TryAdvance(params, data, search.config, search);

  if constexpr (print_progress) {
    static unsigned counter = 0;
    counter++;
    if (counter == 1000000) [[unlikely]] {
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

  if constexpr (debug) std::cout << "Problem: " << problem << std::endl;

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
  std::cout << "Bloom filter population: " << data.bloom->items << std::endl;
  std::cout << "Bloom filter approx    : "
            << data.bloom->ApproximatePopulation() << std::endl;
  std::cout << "Bloom filter error rate: " << data.bloom->ApproximateErrorRate()
            << std::endl;
}
