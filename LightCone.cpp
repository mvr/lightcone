#include <vector>

#include "LifeAPI/LifeAPI.hpp"
#include "LifeAPI/Symmetry.hpp"

#include "Bloom.hpp"
#include "Countdown.hpp"
#include "Params.hpp"

const bool debug = false;
const bool print_progress = true;
const unsigned print_progress_frequency = 1000000;

const bool debug_bloom = false;
const auto debug_bloom_pattern = "";
const auto debug_bloom_key = "";

const unsigned approachRadius = 2; // Needs to match catalyst input
const unsigned placementRequiredLookahead = 4; // This value is borderline unsafe, any larger drops solutions
const unsigned bloomPopulationThreshold = 12; // Min population

const unsigned maxStationaryGens = 64 - 1;

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
  LifeState contact;
  std::vector<Approach> approaches; // TODO: should probably be struct-of-vectors
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
  unsigned ContactRadius() const;
};

CatalystData CatalystData::FromParamsNormal(CatalystParams &params) {
  // LifeHistoryState &approach = params.approaches[0];
  for (auto &approach : params.approaches) {
    approach.AlignWith(params.state);
  }

  LifeState commonContact = ~LifeState();
  for(auto &approach : params.approaches) {
    LifeState reactionWithCatalyst = approach.state;
    reactionWithCatalyst.Step();
    reactionWithCatalyst &= ~params.state;

    LifeState reactionWithoutCatalyst =
        approach.state & approach.marked;
    reactionWithoutCatalyst.Step();

    LifeState contact = reactionWithCatalyst ^ reactionWithoutCatalyst;
    commonContact &= contact;
  }

  std::pair<int, int> contactOrigin;
  unsigned contactCount;
  LifeState tocheck = commonContact;
  for (auto cell = tocheck.FirstOn(); cell != std::make_pair(-1, -1);
       tocheck.Erase(cell), cell = tocheck.FirstOn()) {
    unsigned count =
        (params.approaches[0].state & ~params.state).CountNeighbours(cell);

    bool allequal = std::all_of(
        params.approaches.begin(), params.approaches.end(),
        [&](const auto &approach) {
          unsigned othercount =
              (approach.state & ~params.state).CountNeighbours(cell);
          return othercount == count;
        });

    if (allequal) {
      contactOrigin = cell;
      contactCount = count;
      break;
    }
  }

  params.state.Move(-contactOrigin.first, -contactOrigin.second);
  commonContact.Move(-contactOrigin.first, -contactOrigin.second);
  // TODO: Calculate approach from the soups?
  // This tramples the contents of `params`, seems bad
  for (auto &approach : params.approaches) {
    approach.Move(-contactOrigin.first, -contactOrigin.second);
  }

  params.required.AlignWith(params.state);

  CatalystData result;

  result.contactType =
    contactCount == 1
      ? ContactType::CONTACT1
      : (contactCount == 2 ? ContactType::CONTACT2 : ContactType::CONTACTM);

  result.state = params.state;
  result.halo = params.state.ZOI() & ~params.state;
  result.required = params.required.marked;
  result.state.InteractionCounts(result.history1, result.history2,
                                 result.historyM);
  result.historyFlipped1 = result.history1.Mirrored();
  result.historyFlipped2 = result.history2.Mirrored();
  result.historyFlippedM = result.historyM.Mirrored();
  result.contact = commonContact;
  for (auto &approachPat : params.approaches) {
    Approach approach;
    approach.approachOn = approachPat.marked & approachPat.state;
    approach.approachOff = (approachPat.marked & ~approachPat.state) |
                       result.state | (result.required & ~result.state);
    approach.RecalculateSignature();

    assert((approach.approachOn & approach.approachOff).IsEmpty());
    // More sanity checks?

    result.approaches.push_back(approach);
  }

  for (auto &forbiddenpat : params.forbiddens) {
    Approach approach;
    forbiddenpat.AlignWith(params.state);
    approach.approachOn = forbiddenpat.state & forbiddenpat.marked;
    approach.approachOff = ~forbiddenpat.state & forbiddenpat.marked;
    approach.RecalculateSignature();
    result.forbidden.push_back(approach);
  }

  result.minRecoveryTime = params.minRecoveryTime;
  result.maxRecoveryTime = params.maxRecoveryTime;

  result.transparent = false;

  if constexpr (debug) {
    std::cout << "Loaded catalyst: " << LifeHistoryState(result.state, LifeState(), LifeState::Cell({0, 0})) << std::endl;
    std::cout << "Required: " << LifeHistoryState(result.state, LifeState(), result.required) << std::endl;
    std::cout << "Contact Type: " << result.contactType << std::endl;
    // std::cout << "Params Approach: " << params.approach << std::endl;
    // std::cout << "Params Required: " << params.required << std::endl;
    for (auto &approach : result.approaches) {
      std::cout << "Approach On: " << LifeHistoryState(result.state, LifeState(), approach.approachOn) << std::endl;
      std::cout << "Approach Off: " << LifeHistoryState(result.state, LifeState(),approach.approachOff) << std::endl;
    }
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

  for (auto &forbiddenpat : params.forbiddens) {
    Approach approach;
    forbiddenpat.AlignWith(params.state);
    approach.approachOn = forbiddenpat.state & forbiddenpat.marked;
    approach.approachOff = ~forbiddenpat.state & forbiddenpat.marked;
    approach.RecalculateSignature();
    result.forbidden.push_back(approach);
  }

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

  std::vector<Approach> newApproaches;
  for (auto &oldApproach : approaches) {
    newApproaches.push_back(oldApproach.Transformed(t));
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
    contact.Transformed(t),
    newApproaches,
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

unsigned CatalystData::ContactRadius() const {
  auto [x1,y1,x2,y2] = contact.XYBounds();
  return std::max({std::abs(x1), std::abs(y1), std::abs(x2), std::abs(y1)});
}

struct Placement {
  std::pair<int, int> pos;
  unsigned catalystIx;
  unsigned gen;
};

std::ostream &operator<<(std::ostream &out, const Placement value) {
  return out << "Placing " << value.catalystIx << " at ("
             << value.pos.first << ", " << value.pos.second << ") on gen " << value.gen;
}

// A set of catalyst placements
struct Configuration {
  LifeState state;
  LifeState catalysts;
  LifeState required;

  unsigned numCatalysts;
  unsigned numTransparent;

  unsigned lastInteraction;

  std::vector<Placement> placements;
  std::vector<LifeTarget> targets; // Pre-shifted catalysts
  std::vector<bool> transparent;

  Configuration()
      : state{}, catalysts{}, required{}, numCatalysts{0}, numTransparent{0},
        lastInteraction{0}, placements{}, targets{}, transparent{} {}
};

// Any precomputed data, constant at every node
struct SearchData {
  std::vector<CatalystData> catalysts;
  std::vector<LifeState> collisionMasks;
  unsigned contactRadius;
  LifeBloom *bloom;
  std::vector<Configuration> *allSolutions;
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
  bool allPresent;
  bool nonTransparentPresent;

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
  allPresent = true;
  nonTransparentPresent = true;
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
      if(!config.transparent[i])
        nonTransparentPresent = false;
    }
    i++;
  }

  recoveredTime = allPresent ? (recoveredTime + 1) : 0;
}

std::pair<LifeState, bool> Lookahead::BloomKey(const Configuration &config) const {
  LifeState toHash = state ^ config.catalysts;
  bool valid = nonTransparentPresent
      && (gen > config.lastInteraction + 2)
      // && (config.numCatalysts == 0 || !config.transparent[config.numCatalysts - 1])
      && (toHash.GetPop() > bloomPopulationThreshold);

  return {toHash, valid};
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
    const LifeState &stationaryViolations = stationaryCountdown.finished;
    if(params.maxStationaryCount == 0 || stationaryCountdown.finished.GetPop() > params.maxStationaryCount) {
      // TODO: This is not quite right: this is choosing a violation
      // at random rather than taking the union of them all. But I
      // don't want to store an entire LifeState in the Problem struct
      std::pair<int, int> cell = stationaryViolations.FirstOn();
      if (cell != std::make_pair(-1, -1))
        return {cell, gen, ProblemType::STATIONARY};
    }
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

  if (gen > startTime + params.maxActiveWindowGens && hasInteracted) {
    // TODO: this should be config.catalystsZOI & (config.catalysts ^ state)
    LifeState active = config.catalysts & ~state;
    std::pair<int, int> cell = active.FirstOn();
    if (cell != std::make_pair(-1, -1))
      return {cell, gen, ProblemType::TOO_LONG};
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
      if (seen) {
        if constexpr (debug_bloom) {
          if ((config.state & ~LifeState::ConstantParse(debug_bloom_pattern).Moved(-32,-32)).IsEmpty()) {
            std::cout << "Key here! " << key << std::endl;
          }
        }
        return {{-1, -1}, gen, ProblemType::BLOOM_SEEN};
      }
    }
  }

  return {{-1, -1}, gen, ProblemType::NONE};
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
  case ProblemType::TOO_LONG:
    return LifeState::NZOIAround(cell, gen - currentgen - 1);
  case ProblemType::WINNER:
  case ProblemType::NO_REACTION:
  case ProblemType::BLOOM_SEEN:
    return ~LifeState();
  case ProblemType::NONE:
    __builtin_unreachable();
  }
}

// mvrnote: name?
// Arranged so 0 -> 1 is increasing information
struct CatalystConstraints {
  LifeState tried;            // Permanently ruled out
  LifeState knownUnplaceable; // Can be rescued by earlier placement
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

Problem DetermineProblem(const SearchParams &params, const SearchData &data,
                         const Configuration &config, SearchNode &search, Lookahead &lookahead) {
  while (true) {
    lookahead.Step(config);

    if constexpr (debug) std::cout << "Lookahead to " << lookahead.state << std::endl;

    Problem problem = lookahead.Problem(params, data, config);

    if (problem.type != ProblemType::NONE)
      return problem;

    if (params.useBloomFilter) {
      auto [key, valid] = lookahead.BloomKey(config);

      if constexpr (debug_bloom) {
        if (key == LifeState::ConstantParse(debug_bloom_key).Moved(-32, -32)) {
          std::cout << "Inserted here! " << config.state << std::endl;
        }
      }

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
  if (catalyst.contactType != ContactType::TRANSPARENT && catalyst.contactType != contactType)
    return PlacementValidity::INVALID_CONTACT;

  constexpr LifeState originMask = LifeState::NZOIAround({0, 0}, approachRadius);

  LifeState centered = state.Moved(-p.pos.first, -p.pos.second);

  if (catalyst.approaches.size() > 0) {
    bool anyValid = false;
    bool anyValidAtContact = false;
    for (auto &approach : catalyst.approaches) {
      if (!approach.MatchesSignature(signature))
        continue;

      LifeState mismatches =
          (approach.approachOn & ~centered) | (approach.approachOff & centered);

      if (mismatches.IsEmpty()) {
        anyValid = true;
        anyValidAtContact = true;
      } else {
        if ((originMask & mismatches).IsEmpty()) {
          anyValidAtContact = true;
        } else {
          //
        }
      }
    }

    if (!anyValid && anyValidAtContact)
      return PlacementValidity::FAILED_ELSEWHERE;
    if (!anyValid && !anyValidAtContact)
      return PlacementValidity::FAILED_CONTACT;
  }

  // TODO: should probably check M too...

  // Check whether this catalyst actually would have interacted in a previous
  // generation
  LifeState pastinteractions =
      (catalyst.history1 & historyCount2.Moved(-p.pos.first, -p.pos.second)) |
      (catalyst.history2 & historyCount1.Moved(-p.pos.first, -p.pos.second));
  if (!pastinteractions.IsEmpty()) {
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

  if constexpr (placementRequiredLookahead > 0) {
    LifeState lookahead = centered | catalyst.state;
    lookahead.Step(placementRequiredLookahead);
    if(!(catalyst.required & (lookahead ^ catalyst.state)).IsEmpty())
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

  somePlaceable &= ~(search.history1 & search.history2 & search.historyM);

  for (unsigned gen = search.lookahead.gen; gen < problem.gen; gen++) {
    if constexpr (debug) std::cout << "Gen " << gen << " state: " << current << std::endl;

    if (search.lookahead.hasInteracted && gen > search.lookahead.startTime + params.maxActiveWindowGens)
      break;

    LifeState lightcone = problem.LightCone(gen);
    LifeState lightconeMargin = problem.LightCone(gen - data.contactRadius);
    LifeState possiblePlacements = somePlaceable & lightconeMargin;

    if (possiblePlacements.IsEmpty())
      break;

    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
        currentCountM(UNINITIALIZED);
    current.InteractionCounts(currentCount1, currentCount2, currentCountM);

    LifeState newContactPoints = (currentCount1 & ~currentHistory1) |
                                 (currentCount2 & ~currentHistory2) |
                                 (currentCountM & ~currentHistoryM);

    newContactPoints &= possiblePlacements;

    for (auto cell = newContactPoints.FirstOn(); cell != std::make_pair(-1, -1);
         newContactPoints.Erase(cell), cell = newContactPoints.FirstOn()) {

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
        if (catalyst.contactType != contactType)
          continue;

        auto validSig = std::any_of(catalyst.approaches.begin(),
                                    catalyst.approaches.end(),
                                    [&](const Approach &approach) {
                                      return approach.MatchesSignature(signature);
                                    });
        if (!validSig) {
          search.constraints[i].knownUnplaceable.Set(cell);
          continue;
        }

        if (search.constraints[i].tried.Get(cell) ||
            search.constraints[i].knownUnplaceable.Get(cell)) {
          continue;
        }

        if (!lightcone.Get(cell)) {
          // Need to make sure at least one contact cell is actually in the
          // lightcone
          // TODO: this should be checked by a faster method
          if ((catalyst.contact.Moved(cell) & lightcone).IsEmpty())
            continue;
        }

        Placement p = {cell, i, gen};

        PlacementValidity validity = TestPlacement(
            data, search, current, p, contactType, signature, currentHistory1,
            currentHistory2, currentCount1, currentCount2);

        switch (validity) {
        case PlacementValidity::VALID:
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
      if (search.config.numTransparent < params.maxTransparent && lightcone.Get(cell)) {
        for (unsigned i = 0; i < data.catalysts.size(); i++) {
          const CatalystData &catalyst = data.catalysts[i];
          if (catalyst.contactType != ContactType::TRANSPARENT)
            continue;

          LifeState newContactPoints(UNINITIALIZED);

          switch (contactType /* of the active region */) {
          case ContactType::CONTACT1:
            newContactPoints = catalyst.historyFlipped2.Moved(cell);
            break;
          case ContactType::CONTACT2:
            newContactPoints = catalyst.historyFlipped1.Moved(cell);
            break;
          case ContactType::CONTACTM:
            newContactPoints = (catalyst.historyFlipped1 | catalyst.historyFlipped2).Moved(cell);
            break;
          case ContactType::TRANSPARENT:
            // Not reachable
            break;
          }

          newContactPoints &= ~(search.constraints[i].tried | search.constraints[i].knownUnplaceable);

          for (auto cell = newContactPoints.FirstOn();
               cell != std::make_pair(-1, -1); newContactPoints.Erase(cell),
                    cell = newContactPoints.FirstOn()) {
            Placement p = {cell, i, gen};

            PlacementValidity validity = TestPlacement(
                data, search, current, p, contactType, signature,
                currentHistory1, currentHistory2, currentCount1, currentCount2);

            switch (validity) {
            case PlacementValidity::VALID:
              result.push_back(p);
              break;
            default:
              break;
            }

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
  search.config.transparent.push_back(catalyst.transparent);

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

void RunSearch(const SearchParams &params, const SearchData &data,
               SearchNode &search, Problem problem) {
  if constexpr (debug) std::cout << "Starting node: " << search.config.state << std::endl;

  if constexpr (print_progress) {
    static unsigned counter = 0;
    counter++;
    if (counter == print_progress_frequency) [[unlikely]] {
      std::cout << "Current configuration: " << search.config.state << std::endl;
      LifeState problemGen = search.lookahead.state;
      problemGen.Step(problem.gen - search.lookahead.gen);
      std::cout << "Current problem: " << problem << std::endl;
      // if(problem.cell.first != -1)
      //   std::cout << LifeHistoryState(problemGen, LifeState(), LifeState::Cell(problem.cell)) << std::endl;
      if(params.useBloomFilter) {
        std::cout << "Bloom filter population: " << data.bloom->items << std::endl;
        std::cout << "Bloom filter error rate: " << data.bloom->ApproximateErrorRate() << std::endl;
      }
      // std::cout << "Current placements:" << std::endl;
      // LifeState progression = params.state.state;
      // for (auto &p : search.config.placements) {
      //   const CatalystData &catalystdata = data.catalysts[p.catalystIx];
      //   const LifeState catalyst = catalystdata.state.Moved(p.pos);
      //   progression |= catalyst;
      //   std::cout << p << std::endl;
      //   std::cout << progression << std::endl;
      // }
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

  if (problem.type == ProblemType::TOO_LONG) {
    std::cout << "Too long: " << search.config.state << std::endl;
  }

  if (problem.type == ProblemType::WINNER) {
    std::cout << "Winner: " << search.config.state << std::endl;
    data.allSolutions->push_back(search.config);
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

  std::vector<SearchNode> subsearches;
  subsearches.reserve(placements.size());

  std::vector<Problem> problems;
  problems.reserve(placements.size());

  Lookahead nextPlacementLookahead = search.lookahead;
  LifeState safeContacts1, safeContacts2, safeContactsM;
  bool advanceable = true;

  for (auto &placement : placements) {
    if (search.constraints[placement.catalystIx].tried.Get(placement.pos))
      continue;

    // Placements must be in generation order for this to make sense!
    while (nextPlacementLookahead.gen < placement.gen) {
      LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
        currentCountM(UNINITIALIZED);
      nextPlacementLookahead.state.InteractionCounts(currentCount1, currentCount2, currentCountM);

      if (search.config.numCatalysts > 0) {
        auto lastPlacement = search.config.placements.back();
        if(search.lookahead.gen <= lastPlacement.gen) {
          safeContacts1 |= currentCount1;
          safeContacts2 |= currentCount2;
          safeContactsM |= currentCountM;
        } else {
          LifeState tooClose = LifeState::NZOIAround(placement.pos, approachRadius + search.lookahead.gen - lastPlacement.gen - 1);
          safeContacts1 |= currentCount1 & ~tooClose;
          safeContacts2 |= currentCount2 & ~tooClose;
          safeContactsM |= currentCountM & ~tooClose;
        }
      }

      nextPlacementLookahead.Step(search.config);

      if (advanceable) {
        LifeState missed = (search.lookahead.state ^ search.config.catalysts).ZOI() & ~problem.LightCone(search.lookahead.gen);

        if (!missed.IsEmpty()) {
          advanceable = false;
        } else {
          search.history1 |= currentCount1;
          search.history2 |= currentCount2;
          search.historyM |= currentCountM;
          search.lookahead = nextPlacementLookahead;
          if constexpr (debug) std::cout << "Advanced early to " << search.lookahead.state << std::endl;
        }
      }
    }

    search.constraints[placement.catalystIx].tried.Set(placement.pos);

    subsearches.push_back(search);
    SearchNode &newSearch = subsearches.back();

    MakePlacement(params, data, newSearch, placement);

    // Now invalidate every placement that isn't safe
    // (Doesn't reset the `tried` field)
    for (unsigned i = 0; i < data.catalysts.size(); i++) {
      switch (data.catalysts[i].contactType) {
      case ContactType::CONTACT1:
        newSearch.constraints[i].knownUnplaceable &= safeContacts1;
        break;
      case ContactType::CONTACT2:
        newSearch.constraints[i].knownUnplaceable &= safeContacts2;
        break;
      case ContactType::CONTACTM:
        newSearch.constraints[i].knownUnplaceable &= safeContactsM;
        break;
      case ContactType::TRANSPARENT:
        newSearch.constraints[i].knownUnplaceable = LifeState();
        // They are always placeable
        break;
      }
    }

    Lookahead placed = nextPlacementLookahead;
    placed.state |= data.catalysts[placement.catalystIx].state.Moved(placement.pos);
    placed.missingTime.push_back(0);
    placed.catalystHasInteracted.push_back(false);

    problems.push_back(DetermineProblem(params, data, newSearch.config, newSearch, placed));
  }

  for (unsigned i = 0; i < subsearches.size(); i++) {
    if constexpr (debug) {
      if (params.hasOracle && !(subsearches[i].config.state & ~params.oracle).IsEmpty()) {
        std::cout << "Oracle failed: " << subsearches[i].config.state << std::endl;
        continue;
      }
    }
    if constexpr (debug) std::cout << "Branching node: " << subsearches[i].config.placements.back() << std::endl;

    RunSearch(params, data, subsearches[i], problems[i]);
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

void PrintSummary(std::vector<Configuration> &pats) {
  std::cout << "x = 0, y = 0, rule = B3/S23" << std::endl;
  for (unsigned i = 0; i < pats.size(); i += 8) {
    std::vector<Configuration> rowConfigurations = std::vector<Configuration>(pats.begin() + i, pats.begin() + std::min((unsigned)pats.size(), i + 8));
    std::vector<LifeState> row;
    for (auto &s : rowConfigurations) {
      row.push_back(s.state);
    }
    std::cout << RowRLE(row) << std::endl;
  }
}

int main(int, char *argv[]) {
  auto toml = toml::parse(argv[1]);
  SearchParams params = SearchParams::FromToml(toml);

  std::vector<Configuration> allSolutions;

  std::vector<CatalystData> catalystdata;
  for (auto &c : params.catalysts) {
    auto newdata = CatalystData::FromParams(c);
    catalystdata.insert(catalystdata.end(), newdata.begin(), newdata.end()); // Why is C++ like this
  }

  unsigned contactRadius = 0;
  for (auto &c : catalystdata) {
    contactRadius = std::max(contactRadius, c.ContactRadius());
  }

  std::vector<LifeState> masks = CalculateCollisionMasks(catalystdata);

  LifeBloom *bloom = 0;
  if (params.useBloomFilter)
    bloom = new LifeBloom();

  SearchData data = {catalystdata, masks, contactRadius, bloom, &allSolutions};

  SearchNode search(params, data);

  search.BlockEarlyInteractions(params, data);

  Lookahead lookahead = search.lookahead;
  Problem problem = DetermineProblem(params, data, search.config, search, lookahead);

  RunSearch(params, data, search, problem);

  if (params.printSummary) {
    std::cout << "All solutions:" << std::endl;
    PrintSummary(allSolutions);
  }

  if (params.useBloomFilter) {
    std::cout << "Bloom filter population: " << data.bloom->items << std::endl;
    std::cout << "Bloom filter approx    : "
              << data.bloom->ApproximatePopulation() << std::endl;
    std::cout << "Bloom filter error rate: "
              << data.bloom->ApproximateErrorRate() << std::endl;
  }
}

// Ideas:

// DONE: the main problem is identifying faster when a catalyst
// placement is bad, without having to re-run the lookahead all the
// way up to that placement point. Is it correct that when we make a
// placement, the new problem is *always* after that placement is
// reached? Then we can pass the state lookahead down to the child
// edit: this helped a little but not much

// TODO: do a kind of fast pass in `DetermineProblem` for `REQUIRED`
// violations. this will mean needing to calculate the earliest that a
// non-`REQUIRED` problem that can happen, and restarting if we reach
// that point

// TODO: maybe during collect placements, run a life-with-unknowns
// (together with / instead of influencable) and see whether anything
// *past* the problem is inevitable

// This seems like it could happen at branching time rather than
// collection time, we can keep track of "influencable" cells and do a
// short lookahead at each placement to see if `required` is broken.


// TODO: use a perfect hash function for the signatures, so they can
// all be checked quickly (checking signatures may not have a time cost worth
// worrying about)

// TODO: I have probably killed performance with a lot of these
// changes... need to do some profiling

// TODO: It might be possible to check the `Contains` in
// `Lookahead::Step` more efficiently. Currently it does a full-board
// calculation for each placed catalyst, but most catalysts spend most
// of their time recovered. How about it instead looks at the active
// cells within the ZOI of a catalyst, and then iterates through those
// active cells and sees which catalyst they correspond to.

// TODO: debug print when there is a chain of placements that are in
// reverse generation order. How often does that kind of cascade
// actually happen?

// TODO: how about checking to see whether the earlier problem still
// exists, and re-using that problem if so?

// TODO: should transparent catalysts be split into many copies,
// differing in which cell gets hit first? This might streamline some
// of their handling

// DONE: argh, I really wanted to test just one contact cell per
// catalyst. But we can get unlucky and choose the wrong contact cell,
// and interactions right on the edge of the light cone can get
// dropped.

// DONE: ResetLightcone could be merged into the main loop, it may be
// worth it

// TODO: there are some very small LifeStates that could be stored/queried
// more efficiently by just storing the cell coordinates
