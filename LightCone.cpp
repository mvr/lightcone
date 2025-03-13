#include <iostream>
#include <vector>

#include "LifeAPI/LifeAPI.hpp"
#include "LifeAPI/LifeTarget.hpp"
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

  static Approach FromParams(const LifeHistory &params);
};

Approach Approach::FromParams(const LifeHistory &params) {
  LifeState catalyst = params.state & ~params.marked;

  // LifeState reactionWithCatalyst = params.state.Stepped() & ~catalyst;
  // LifeState reactionWithoutCatalyst = (params.state & params.marked).Stepped();

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
  bool limited;

  std::vector<Approach> forbidden;

  static CatalystData FromParamsNormal(CatalystParams &params);
  static CatalystData FromParamsTransparent(CatalystParams &params);
  static std::vector<CatalystData> FromParams(CatalystParams &params);

  CatalystData Transformed(SymmetryTransform t);

  LifeState CollisionMask(const CatalystData &b) const;
  unsigned ContactRadius() const;
};

CatalystData CatalystData::FromParamsNormal(CatalystParams &params) {
  for (auto &approach : params.approaches) {
    approach.AlignWith(params.state);
  }

  LifeState commonContact = ~LifeState();
  for(auto &approach : params.approaches) {
    LifeState reactionWithCatalyst = approach.state.Stepped();
    reactionWithCatalyst &= ~params.state;

    LifeState reactionWithoutCatalyst =
      (approach.state & approach.marked).Stepped();

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
  for (auto &soup : params.soups) {
    soup.Move(-contactOrigin.first, -contactOrigin.second);
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
  result.limited = params.limited;

  if constexpr (debug) {
    std::cout << "Loaded catalyst: " << LifeHistory(result.state, LifeState(), LifeState::Cell({0, 0})) << std::endl;
    std::cout << "Required: " << LifeHistory(result.state, LifeState(), result.required) << std::endl;
    std::cout << "Contact Type: " << result.contactType << std::endl;
    // std::cout << "Params Approach: " << params.approach << std::endl;
    // std::cout << "Params Required: " << params.required << std::endl;
    for (auto &approach : result.approaches) {
      std::cout << "Approach On: " << LifeHistory(result.state, LifeState(), approach.approachOn) << std::endl;
      std::cout << "Approach Off: " << LifeHistory(result.state, LifeState(),approach.approachOff) << std::endl;
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
  result.limited = params.limited;

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
    limited,
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

  // Block active intersecting with required
  result |= (state.ZOI() & ~required).Convolve(b.required.Mirrored());

  // And vice versa
  result |= (b.state.ZOI() & ~b.required).Mirrored().Convolve(required);

  LifeState newresult;

  LifeState selfActive = (state.ZOI() & ~required) | state;
  // LifeState activeCount1(UNINITIALIZED), activeCount2(UNINITIALIZED), activeCountM(UNINITIALIZED);
  // selfActive.InteractionCounts(activeCount1, activeCount2, activeCountM);

  // newresult |= activeCount2.Convolve((b.required & b.history1).Mirrored());
  // newresult |= activeCount1.Convolve((b.required & b.history2).Mirrored());

  newresult |= selfActive.Convolve(b.contact.Mirrored());

  LifeState bActive = (b.state.ZOI() & ~b.required) | b.state;
  // LifeState bActiveCount1(UNINITIALIZED), bActiveCount2(UNINITIALIZED), bActiveCountM(UNINITIALIZED);
  // bActive.InteractionCounts(bActiveCount1, bActiveCount2, bActiveCountM);

  // newresult |= bActiveCount2.Mirrored().Convolve(required & history1);
  // newresult |= bActiveCount1.Mirrored().Convolve(required & history2);

  newresult |= bActive.Mirrored().Convolve(contact);

  return result | newresult;
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
  unsigned numLimited;

  unsigned lastInteraction;

  std::vector<Placement> placements;
  std::vector<LifeTarget> targets; // Pre-shifted catalysts
  std::vector<bool> transparent;

  Configuration()
      : state{}, catalysts{}, required{}, numCatalysts{0}, numTransparent{0}, numLimited{0},
        lastInteraction{0}, placements{}, targets{}, transparent{} {}

  void MakePlacement(const SearchParams &params, const CatalystData &catalyst,
                     const Placement &placement);
};

void Configuration::MakePlacement(const SearchParams &params, const CatalystData &catalyst,
                   const Placement &placement) {
  const LifeState catalystState = catalyst.state.Moved(placement.pos);

  numCatalysts++;
  if(catalyst.transparent) numTransparent++;
  if(catalyst.limited) numLimited++;
  lastInteraction = std::max(lastInteraction, placement.gen);
  state |= catalystState;
  catalysts |= catalystState;
  required |= catalyst.required.Moved(placement.pos);
  placements.push_back(placement);
  targets.push_back(
      LifeTarget(catalyst.state.Moved(placement.pos),
                 catalyst.halo.Moved(placement.pos)));
  transparent.push_back(catalyst.transparent);
}

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
  REQUIRED,
  FILTER,
  UNRECOVERED,
  TOO_LONG,
  NO_REACTION,
  NOT_TRANSPARENT,
  STATIONARY,
};

std::ostream &operator<<(std::ostream &out, const ProblemType value) {
  return out << [value]() {
    switch (value) {
    case ProblemType::NONE:            return "NONE";
    case ProblemType::REQUIRED:        return "REQUIRED";
    case ProblemType::FILTER:          return "FILTER";
    case ProblemType::UNRECOVERED:     return "UNRECOVERED";
    case ProblemType::TOO_LONG:        return "TOO_LONG";
    case ProblemType::NO_REACTION:     return "NO_REACTION";
    case ProblemType::NOT_TRANSPARENT: return "NOT_TRANSPARENT";
    case ProblemType::STATIONARY:      return "STATIONARY";
    }
  }();
}

struct Problem {
  std::pair<int, int> cell;
  unsigned gen;
  ProblemType type;

  LifeState LightCone(unsigned gen) const;
  bool IsGlobal() const;
  bool Subsumes(Problem other) const;
};

std::ostream &operator<<(std::ostream &out, const Problem value) {
  return out << value.type << " on gen " << value.gen << " at ("
             << value.cell.first << ", " << value.cell.second << ")";
}

struct LookaheadOutcome {
  Problem problem;
  unsigned bloomSeenGen;
  unsigned timeoutGen;
  unsigned winnerGen;
  bool winner;
};


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

  void MakePlacement(const SearchParams &params, const CatalystData &catalyst,
                     const Placement &placement);
  
  void Step(const Configuration &config);
  Problem CurrentProblem(const SearchParams &params, const SearchData &data,
                         const Configuration &config) const;

  std::pair<LifeState, bool> BloomKey(const Configuration &config) const;
};

void Lookahead::MakePlacement(const SearchParams &params, const CatalystData &catalyst,
                              const Placement &placement) {
  state |= catalyst.state.Moved(placement.pos);
  missingTime.push_back(0);
  catalystHasInteracted.push_back(false);
}

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

Problem Lookahead::CurrentProblem(const SearchParams &params, const SearchData &data,
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

  return {{-1, -1}, gen, ProblemType::NONE};
}

// Contact points that are close enough to the current problem to
// have an effect on it.
LifeState Problem::LightCone(unsigned currentgen) const {
  switch (type) {
  case ProblemType::REQUIRED:
  case ProblemType::FILTER:
  case ProblemType::UNRECOVERED:
  case ProblemType::NOT_TRANSPARENT:
  case ProblemType::STATIONARY:
  case ProblemType::TOO_LONG:
    return LifeState::NZOIAround(cell, std::abs((int)gen - (int)currentgen - 1));
  case ProblemType::NO_REACTION:
  case ProblemType::NONE:
    return ~LifeState();
  }
}

bool Problem::IsGlobal() const {
  switch (type) {
  case ProblemType::REQUIRED:
  case ProblemType::FILTER:
  case ProblemType::UNRECOVERED:
  case ProblemType::NOT_TRANSPARENT:
  case ProblemType::STATIONARY:
  case ProblemType::TOO_LONG:
    return false;
  case ProblemType::NO_REACTION:
  case ProblemType::NONE:
    return true;
  }
}


// Sooner and smaller
bool Problem::Subsumes(Problem other) const {
  if (gen >= other.gen) return false;

  if (IsGlobal())
    return false;

  if (!IsGlobal() && other.IsGlobal())
    return true;

  // Now neither is global
  return other.LightCone(gen).Get(cell);
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

LookaheadOutcome DetermineProblem(const SearchParams &params, const SearchData &data,
                                  const Configuration &config, Lookahead &lookahead) {
  unsigned bloomSeenGen = std::numeric_limits<unsigned>::max();
  unsigned timeoutGen = std::numeric_limits<unsigned>::max();
  unsigned winnerGen = std::numeric_limits<unsigned>::max();
  bool winner = false;

  while (true) {
    lookahead.Step(config);

    if constexpr (debug) std::cout << "Lookahead to " << lookahead.state << std::endl;

    Problem problem = lookahead.CurrentProblem(params, data, config);

    if (!winner &&
        lookahead.gen > params.minFirstActiveGen && lookahead.hasInteracted &&
        lookahead.recoveredTime >= params.minStableTime) {
      winnerGen = lookahead.gen;
      winner = true;
    }

    if (lookahead.hasInteracted && lookahead.recoveredTime > params.maxStableTime)
      timeoutGen = lookahead.gen;

    bool shouldReturn = problem.type != ProblemType::NONE ||
                        (!params.continueAfterSuccess && winner) ||
                        (lookahead.hasInteracted && lookahead.recoveredTime > params.maxStableTime);

    if (shouldReturn)
      return LookaheadOutcome(problem, bloomSeenGen, timeoutGen, winnerGen, winner);

    if (params.useBloomFilter && bloomSeenGen == std::numeric_limits<unsigned>::max()) {
      auto [key, valid] = lookahead.BloomKey(config);

      if(valid) {
        bool seen = data.bloom->Lookup(key);
        if (seen) {
          if constexpr (debug_bloom) {
            if ((config.state & ~LifeState::ConstantParse(debug_bloom_pattern).Moved(-32,-32)).IsEmpty()) {
              std::cout << "Key here! " << key << std::endl;
            }
          }

          bloomSeenGen = lookahead.gen;
        } else {
          if constexpr (debug_bloom) {
            if (key == LifeState::ConstantParse(debug_bloom_key).Moved(-32, -32)) {
              std::cout << "Inserted here! " << config.state << std::endl;
            }
          }

          data.bloom->Insert(key);
        }
      }
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

  if (catalyst.approaches.size() > 0) [[likely]] {
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

    if (!anyValid && !anyValidAtContact)
      return PlacementValidity::FAILED_CONTACT;
    if (!anyValid && anyValidAtContact)
      return PlacementValidity::FAILED_ELSEWHERE;
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
  if (!immediatebirths.IsEmpty()) {
    if ((originMask & immediatebirths).IsEmpty()) {
      return PlacementValidity::FAILED_ELSEWHERE;
    } else {
      return PlacementValidity::FAILED_CONTACT;
    }
  }

  // Check the forbidden approaches
  for (auto &f : catalyst.forbidden) {
    LifeState differences = (f.approachOn & ~centered) |
                            (f.approachOff & centered);
    if (differences.IsEmpty())
      return PlacementValidity::FAILED_ELSEWHERE;
  }

  // Check whether the catalyst explodes very quickly
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
                                         SearchNode &search, LookaheadOutcome &constraint) {
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

  // TODO: this will exclude some transparent placements, no?
  LifeState somePlaceable = LifeState();
  for (unsigned i = 0; i < data.catalysts.size(); i++) {
    somePlaceable |= ~(search.constraints[i].tried | search.constraints[i].knownUnplaceable);
  }

  somePlaceable &= ~(search.history1 & search.history2 & search.historyM);

  unsigned lastPlacementGen = std::min(std::min(constraint.bloomSeenGen, constraint.timeoutGen), constraint.problem.gen);

  for (unsigned gen = search.lookahead.gen; gen < lastPlacementGen; gen++) {
    if constexpr (debug) std::cout << "Gen " << gen << " state: " << current << std::endl;

    if (search.lookahead.hasInteracted && gen > search.lookahead.startTime + params.maxActiveWindowGens)
      break;

    LifeState lightcone = constraint.problem.LightCone(gen);
    LifeState lightconeMargin = constraint.problem.LightCone(gen - data.contactRadius);
    LifeState possiblePlacements = somePlaceable & lightconeMargin;

    if (possiblePlacements.IsEmpty())
      break;

    LifeState currentCount1(InitializedTag::UNINITIALIZED), currentCount2(InitializedTag::UNINITIALIZED),
        currentCountM(InitializedTag::UNINITIALIZED);
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

        if (catalyst.limited && search.config.numLimited == params.maxLimited)
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

          LifeState newContactPoints(InitializedTag::UNINITIALIZED);

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

  search.config.MakePlacement(params, catalyst, placement);
  search.lookahead.MakePlacement(params, catalyst, placement);

  search.history1 |= catalyst.history1.Moved(placement.pos);
  search.history2 |= catalyst.history2.Moved(placement.pos);
  search.historyM |= catalyst.historyM.Moved(placement.pos);

  for (unsigned t = 0; t < data.catalysts.size(); t++) {
    if(data.catalysts[t].limited && search.config.numLimited == params.maxLimited) continue;
    search.constraints[t].tried |=
        data.collisionMasks[placement.catalystIx * data.catalysts.size() + t]
            .Moved(placement.pos.first, placement.pos.second);
  }
}

void RunSearch(const SearchParams &params, const SearchData &data,
               SearchNode &search, LookaheadOutcome constraint) {
  if constexpr (debug) std::cout << "Starting node: " << search.config.state << std::endl;

  if constexpr (print_progress) {
    static unsigned counter = 0;
    counter++;
    if (counter == print_progress_frequency) [[unlikely]] {
      std::cout << "Current configuration: " << search.config.state << std::endl;
      LifeState problemGen = search.lookahead.state;
      problemGen.Step(constraint.problem.gen - search.lookahead.gen);
      std::cout << "Current problem: " << constraint.problem << std::endl;
      // if(problem.cell.first != -1)
      //   std::cout << LifeHistory(problemGen, LifeState(), LifeState::Cell(problem.cell)) << std::endl;
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
    std::cout << "Problem: " << constraint.problem << std::endl;
    if (constraint.problem.cell.first != -1) {
      LifeState problemGen = search.lookahead.state;
      problemGen.Step(constraint.problem.gen - search.lookahead.gen);
      std::cout << LifeHistory(problemGen, LifeState(), LifeState::Cell(constraint.problem.cell)) << std::endl;
    }
  }

  if (constraint.problem.type == ProblemType::TOO_LONG) {
    std::cout << "Too long: " << search.config.state << std::endl;
  }

  if (constraint.winner && constraint.winnerGen <= constraint.bloomSeenGen) {
    std::cout << "Winner: " << search.config.state << std::endl;
    data.allSolutions->push_back(search.config);
    if constexpr (debug) {
      std::cout << "Required: " << LifeHistory(search.config.state, LifeState(), search.config.required) << std::endl;
      for (auto &p : search.config.placements) {
        std::cout << "Placement: " << p.catalystIx << " at (" << p.pos.first << ", " << p.pos.second << ")" << std::endl;
      }
    }
    if (!params.continueAfterSuccess)
      return;
  }

  if (search.config.numCatalysts == params.maxCatalysts)
    return;

  std::vector<Placement> placements =
      CollectPlacements(params, data, search, constraint);

  std::vector<LookaheadOutcome> problems;
  problems.reserve(placements.size());

  Lookahead nextPlacementLookahead = search.lookahead;

  for (const auto &placement : placements) {
    // Placements in the list must be in generation order for this to make sense!
    while (nextPlacementLookahead.gen < placement.gen) {
      LifeState currentCount1(InitializedTag::UNINITIALIZED), currentCount2(InitializedTag::UNINITIALIZED),
        currentCountM(InitializedTag::UNINITIALIZED);
      nextPlacementLookahead.state.InteractionCounts(currentCount1, currentCount2, currentCountM);
      nextPlacementLookahead.Step(search.config);
    }

    const CatalystData &catalyst = data.catalysts[placement.catalystIx];
    
    Configuration newConfig = search.config;
    newConfig.MakePlacement(params, catalyst, placement);

    Lookahead newLookahead = nextPlacementLookahead;
    newLookahead.MakePlacement(params, catalyst, placement);

    problems.push_back(DetermineProblem(params, data, newConfig, newLookahead));
  }

  bool advanceable = true;
  nextPlacementLookahead = search.lookahead;  
  LifeState safeContacts1, safeContacts2, safeContactsM;

  unsigned placementIx = 0;
  for (const auto &placement : placements) {
    if (search.constraints[placement.catalystIx].tried.Get(placement.pos)) {
      placementIx++;
      continue;
    }
    search.constraints[placement.catalystIx].tried.Set(placement.pos);
    
    // Placements in the list must be in generation order for this to make sense!
    while (nextPlacementLookahead.gen < placement.gen) {
      LifeState currentCount1(InitializedTag::UNINITIALIZED), currentCount2(InitializedTag::UNINITIALIZED),
          currentCountM(InitializedTag::UNINITIALIZED);
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
        LifeState missed = (search.lookahead.state ^ search.config.catalysts).ZOI() & ~constraint.problem.LightCone(search.lookahead.gen);

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

    SearchNode newSearch = search;

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

    if constexpr (debug) {
      if (params.hasOracle && !(newSearch.config.state & ~params.oracle).IsEmpty()) {
        std::cout << "Oracle failed: " << newSearch.config.state << std::endl;
        placementIx++;
        continue;
      }
    }
    if constexpr (debug) std::cout << "Branching node: " << newSearch.config.placements.back() << std::endl;

    RunSearch(params, data, newSearch, problems[placementIx]);

    placementIx++;
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
    LifeState currentCount1(InitializedTag::UNINITIALIZED), currentCount2(InitializedTag::UNINITIALIZED),
        currentCountM(InitializedTag::UNINITIALIZED);
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

  if (params.maxStationaryTime != 0 &&
      (unsigned)params.maxStationaryTime > maxStationaryGens) {
    std::cout << "`max-stationary-time` is higher than allowed by the hardcoded value!" << std::endl;
    exit(1);
  }

  std::vector<Configuration> allSolutions;

  std::vector<CatalystData> catalystdata;
  for (auto &c : params.catalysts) {
    auto newdata = CatalystData::FromParams(c);
    catalystdata.insert(catalystdata.end(), newdata.begin(), newdata.end()); // Why is C++ like this
  }

  unsigned contactRadius = 0;
  for (auto &c : catalystdata) {
    if(c.transparent) continue;
    contactRadius = std::max(contactRadius, c.ContactRadius());
  }
  if constexpr (debug) std::cout << "Contact radius: " << contactRadius << std::endl;

  std::vector<LifeState> masks = CalculateCollisionMasks(catalystdata);

  LifeBloom *bloom = 0;
  if (params.useBloomFilter)
    bloom = new LifeBloom();

  SearchData data = {catalystdata, masks, contactRadius, bloom, &allSolutions};

  SearchNode search(params, data);

  search.BlockEarlyInteractions(params, data);

  Lookahead lookahead = search.lookahead;
  LookaheadOutcome problem = DetermineProblem(params, data, search.config, lookahead);

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
// *past* the problem is inevitable. This seems like it could happen
// at branching time rather than collection time, we can keep track of
// "influencable" cells and do a short lookahead at each placement to
// see if `required` is broken. On reflection, the lightcone placement
// reasoning should basically handle this

// TODO: use a perfect hash function for the signatures, so they can
// all be checked quickly (checking signatures may not have a time
// cost worth worrying about)

// TODO: I have probably killed performance with a lot of these
// changes... need to do some profiling

// TODO: Check the `Contains` in `Lookahead::Step` more efficiently.
// Currently it does a full-board calculation for each placed
// catalyst, but most catalysts spend most of their time recovered.
// How about it instead looks at the active cells within the ZOI of a
// catalyst, and then iterates through those active cells and sees
// which catalyst they correspond to.

// TODO: Identify when there is a chain of placements that are in
// reverse generation order. How often does that kind of cascade
// actually happen?

// TODO: Should transparent catalysts be split into many copies,
// differing in which cell gets hit first? This might streamline some
// of their handling, but we would have to be more careful about
// avoiding identical placements coming from different copies.

// DONE: I really wanted to test just one contact cell per catalyst.
// But we can get unlucky and choose the wrong contact cell, and
// interactions right on the edge of the lightcone can get dropped.

// DONE: `ResetLightcone` could be merged into the main loop, it may be
// worth it

// TODO: There are some very small `LifeState`s that could be
// stored/queried more efficiently by just storing a list of cell
// coordinates

// TODO: Allow non-transparent catalysts that don't supply an approach

// TODO: How about checking to see whether the earlier problem still
// exists, and re-using that problem if so? There should be a
// threshold where if a new placement explodes fast enough, we use
// that new problem, otherwise we use the old one. Otherwise, we place
// a bunch of hopeless catalysts to solve the original problem. Or, if
// the new placement causes a problem contained in the lightcone of
// the previous one, then there is no harm solving that more recent
// problem first. After trying this, there's an issue: it will place a
// dense field of failing catalysts, all contained within the
// lightcone of the original problem. We will need some strategy like
// avoiding placements that are within the lightcone of known
// problems. Or not considering a problem pressing if it swamped by an
// earlier problem.

// TODO: Another idea is to identify the problem with the fewest
// remaining interaction points in its lightcone, and then try to
// solve that one. Or could it be good enough to just identify when an
// early problem has had all its interaction points run out? After we
// have all the placements, we can identify whether the problem spots
// are still influencable

// TODO: Or maybe, we prefer problems that occur in the lightcone of
// the previous problem.

// DONE: We could add a `max-stationary-count` constraint that allows
// some number of cells to be stationary, rather than forbidding them
// all

// TODO: Reordering placements to prefer nearby ones could help with
// the quadratic explosion with large active region. Or not: it
// doesn't help starting with later placement first, because those
// placements may not exist when you place something earlier. But what
// if we leverage `ResetPlacements` some more, and have later
// placements knownUnplaceable until a lightcone reaches them? Seems
// complicated

// Another downside is not being able to advance the lookahead as
// often. It might also help the bloom filter out?

// TODO: It would be good if glider-eating reactions could be
// forbidden, but only after their first occurrence along a ray. (So we
// can still eat the glider, but not duplicate work with more distance
// placements.)

// DONE: We should have a `continue-after-success` option, I am
// worried we are missing useful solutions by needing a high
// `min-stable-time`. This would mean both reporting winner status and
// the next problem in `DetermineProblem`

// TODO: What if each subtree returns some information about the
// problems that were encountered. Could that be useful?
