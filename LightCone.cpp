#include <vector>

#include "LifeAPI/LifeAPI.hpp"
#include "LifeAPI/Symmetry.hpp"

#include "Params.hpp"

const bool debug = false;
const bool print_progress = true;

const unsigned approachRadius = 1; // Needs to match catalyst input
// const unsigned perturbationLookahead = 5; // TODO

enum ContactType {
  CONTACT1,
  CONTACT2,
  CONTACTM,
};

std::ostream& operator<<(std::ostream& out, const ContactType value){
  return out << [value]() {
    switch (value) {
    case ContactType::CONTACT1: return "CONTACT1";
    case ContactType::CONTACT2: return "CONTACT2";
    case ContactType::CONTACTM: return "CONTACTM";
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
  ContactType contactType; // The at the origin contact cell
  unsigned maxRecoveryTime;
  unsigned minRecoveryTime;

  static CatalystData FromParams(CatalystParams &params);
  static std::vector<CatalystData> FromParamsSymmetrically(CatalystParams &params);

  CatalystData Transformed(SymmetryTransform t);

  LifeState CollisionMask(const CatalystData &b) const;
};

CatalystData CatalystData::FromParams(CatalystParams &params) {
  LifeState contact;
  {
    LifeState withCatalyst = params.approach.state;
    withCatalyst.Step();

    LifeState withoutCatalyst = params.approach.state & params.approach.marked;
    withoutCatalyst.Step();

    contact = withCatalyst ^ withoutCatalyst;
  }

  auto contactorigin = contact.FirstOn();
  // TODO: Calculate approach from the soups?
  params.approach.Move(-contactorigin.first, -contactorigin.second);

  params.state.AlignWith(params.approach.state & ~params.approach.marked);
  params.required.AlignWith(params.state);

  CatalystData result;

  result.state = params.state;
  result.halo = params.state.ZOI() & ~params.state;
  result.required = params.required.marked;
  result.state.InteractionCounts(result.history1, result.history2,
                                 result.historyM);
  result.approachOn = params.approach.marked & params.approach.state;
  result.approachOff = (params.approach.marked & ~params.approach.state) | result.state;

  unsigned contactCount = result.approachOn.CountNeighbours({0, 0});
  result.contactType =
      contactCount == 1
          ? ContactType::CONTACT1
          : (contactCount == 2 ? ContactType::CONTACT2 : ContactType::CONTACTM);

  result.minRecoveryTime = params.minRecoveryTime;
  result.maxRecoveryTime = params.maxRecoveryTime;

  if constexpr (debug) {
    std::cout << "Loaded catalyst: " << result.state << std::endl;
    // std::cout << "Params Approach: " << params.approach << std::endl;
    // std::cout << "Params Required: " << params.required << std::endl;
    std::cout << "Approach: " << LifeHistoryState(result.state | result.approachOn, LifeState(), result.approachOn | result.approachOff) << std::endl;
    std::cout << "Required: " << LifeHistoryState(result.state, LifeState(), result.required) << std::endl;
    std::cout << "Contact Type: " << result.contactType << std::endl;
  }

  return result;
}

CatalystData CatalystData::Transformed(SymmetryTransform t) {
  return {state.Transformed(t),
          halo.Transformed(t),
          required.Transformed(t),
          history1.Transformed(t),
          history2.Transformed(t),
          historyM.Transformed(t),
          approachOn.Transformed(t),
          approachOff.Transformed(t),
          contactType,
          maxRecoveryTime,
          minRecoveryTime};
}

std::vector<CatalystData> CatalystData::FromParamsSymmetrically(CatalystParams &params) {
  CatalystData data = CatalystData::FromParams(params);
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
};


struct Placement {
  std::pair<int, int> pos;
  unsigned catalystIx;
};

struct Perturbation {
  uint64_t hash;
  std::vector<Placement> placements;
};

// A set of catalyst placements
// mvrnote: name? Solution?
struct Configuration {
  LifeState state;
  LifeState catalysts;
  LifeState required;

  unsigned numCatalysts;

  std::vector<Placement> placements;
  std::vector<LifeTarget> targets; // Pre-shifted catalysts

  Configuration()
      : state{}, catalysts{}, required{}, numCatalysts{0},
        placements{}, targets{} {}
};

// The state of a configuration after stepping
struct Lookahead {
  LifeState state;

  unsigned gen;
  bool hasInteracted;
  std::vector<unsigned> missingTime;
  unsigned recoveredTime;

  Lookahead()
      : state{}, gen{0}, hasInteracted{false},
        missingTime{}, recoveredTime{0} {}

  void Step(const Configuration &config);
};

void Lookahead::Step(const Configuration &config) {
  state.Step();
  gen++;

  // Why is there no easy way to iterate with index? My kingdom for a `zip`
  unsigned i = 0;
  bool allPresent = true;
  for (auto &t : config.targets) {
    if (state.Contains(t)) {
      missingTime[i] = 0;
    } else {
      missingTime[i]++;
      hasInteracted = true;
      allPresent = false;
    }
    i++;
  }

  recoveredTime = allPresent ? (recoveredTime + 1) : 0;
}

enum struct ProblemType {
  NONE,
  REQUIRED,
  FILTER,
  UNRECOVERED,
  NO_REACTION,
};

std::ostream& operator<<(std::ostream& out, const ProblemType value){
  return out << [value]() {
    switch (value) {
    case ProblemType::NONE:        return "NONE";
    case ProblemType::REQUIRED:    return "REQUIRED";
    case ProblemType::FILTER:      return "FILTER";
    case ProblemType::UNRECOVERED: return "UNRECOVERED";
    case ProblemType::NO_REACTION: return "NO_REACTION";
    }
  }();
}

struct Problem {
  std::pair<int, int> cell;
  unsigned gen;
  ProblemType type;

  LifeState LightCone(unsigned gen);
};

std::ostream& operator<<(std::ostream& out, const Problem value){
  return out << value.type << " on gen " << value.gen << " at (" << value.cell.first << ", " << value.cell.second << ")";
}

// Contact points that are close enough to the current problem to
// have an effect on it.
LifeState Problem::LightCone(unsigned currentgen) {
  switch (type) {
  case ProblemType::REQUIRED:
  case ProblemType::FILTER:
  case ProblemType::UNRECOVERED:
    return LifeState::NZOIAround(cell, gen - currentgen);
  case ProblemType::NO_REACTION:
    return ~LifeState();
  case ProblemType::NONE:
    __builtin_unreachable();
  }
}

Problem DetermineProblem(const SearchParams &params, const SearchData &data,
                         const Configuration &config, const Lookahead &start) {
  Lookahead lookahead = start;

  while (true) {
    {
      LifeState requiredViolations =
          config.required &
          (lookahead.state ^ config.catalysts);
      std::pair<int, int> cell = requiredViolations.FirstOn();
      if (cell != std::make_pair(-1, -1))
        return {cell, lookahead.gen, ProblemType::REQUIRED};
    }

    {
      for (unsigned i = 0; i < config.numCatalysts; i++) {
        if (lookahead.missingTime[i] > data.catalysts[config.placements[i].catalystIx].maxRecoveryTime) {
          const LifeTarget &target = config.targets[i];
          std::pair<int, int> cell =
              (target.wanted & ~lookahead.state).FirstOn();
          if (cell.first == -1 && cell.second == -1)
            cell = (target.unwanted & lookahead.state).FirstOn();
          return {cell, lookahead.gen, ProblemType::UNRECOVERED};
        }
      }
    }

    {
      if (lookahead.gen > params.maxFirstActiveGen && !lookahead.hasInteracted)
        return {{-1, -1}, lookahead.gen, ProblemType::NO_REACTION};
    }

    {
      if (lookahead.gen > params.minFirstActiveGen && lookahead.hasInteracted && lookahead.recoveredTime > params.minStableTime) {
        return {{-1, -1}, 0, ProblemType::NONE};
      }
    }

    lookahead.Step(config);
  }
}

// mvrnote: name?
// Arranged so 0 -> 1 is increasing information
struct CatalystConstraints {
  LifeState tried;
  // These refer to the validity *at the origin of the catalyst*
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
  // std::vector<Problem> problems;

  SearchNode(const SearchParams &params, const SearchData &data) {
    config = Configuration();
    lookahead = Lookahead();

    config.state = params.state.state;
    lookahead.state = params.state.state;

    config.state.InteractionCounts(history1, history2,
                                             historyM);
    constraints = std::vector<CatalystConstraints>(data.catalysts.size(),
                                                   CatalystConstraints());
    // TODO: restrict locations catalysts can be placed
  }

  void Step(const SearchParams &params, const SearchData &data);
};

struct ContactEnvelope {
  LifeState firstActive1;
  LifeState firstActive2;
  LifeState firstActiveM;
};

// ContactEnvelope CalculateContactEnvelope(SearchParams &params, SearchNode &state) {
// }

// void SearchNode::Step(const SearchParams &params, const SearchData &data) {
// }

// Wait, do we want to do this? A cell could e.g. be history1 and then history2 later
void BlockIncorrectContacts(const SearchParams &params, const SearchData &data,
                            SearchNode &search, LifeState state) {
  LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
      currentCountM(UNINITIALIZED);
  state.InteractionCounts(currentCount1, currentCount2,
                                      currentCountM);
  for (unsigned i = 0; i < data.catalysts.size(); i++) {
    switch (data.catalysts[i].contactType) {
    case CONTACT1:
      search.constraints[i].tried |= currentCount2 | currentCountM;
      break;
    case CONTACT2:
      search.constraints[i].tried |= currentCount1 | currentCountM;
      break;
    case CONTACTM:
      search.constraints[i].tried |= currentCount1 | currentCount2;
      break;
    }
  }
}

// TODO: determine exactly when this needs to be run.
// Is it after the final perturbation?
void TryAdvance(const SearchParams &params, const SearchData &data,
                SearchNode &search) {
  while (true) {
    if constexpr (debug) std::cout << "Trying to advance: " << search.lookahead.state << std::endl;

    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
        currentCountM(UNINITIALIZED);
    search.lookahead.state.InteractionCounts(
        currentCount1, currentCount2, currentCountM);
    LifeState newContactPoints = (currentCount1 & ~search.history1) |
                                 (currentCount2 & ~search.history2) |
                                 (currentCountM & ~search.historyM);

    if constexpr (debug) std::cout << "New contact points: " << newContactPoints << std::endl;

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
      if constexpr (debug) std::cout << "Advancing!" << std::endl;
      search.lookahead.Step(search.config);
      search.history1 |= currentCount1;
      search.history2 |= currentCount2;
      search.historyM |= currentCountM;
      // BlockIncorrectContacts(params, data, search, search.lookahead.state);
    } else {
      break;
    }
  }
}

enum struct PlacementValidity {
  VALID,
  INVALID_CONTACT,
  FAILED_CONTACT,  // Invalid within `approachRadius` of the origin
  FAILED_ELSEWHERE // Invalid outside that, so will need to be re-checked
};

std::ostream& operator<<(std::ostream& out, const PlacementValidity value){
  return out << [value]() {
    switch (value) {
    case PlacementValidity::VALID: return "VALID";
    case PlacementValidity::INVALID_CONTACT: return "INVALID_CONTACT";
    case PlacementValidity::FAILED_CONTACT: return "FAILED_CONTACT";
    case PlacementValidity::FAILED_ELSEWHERE: return "FAILED_ELSEWHERE";
    }
  }();
}

PlacementValidity TestPlacement(const SearchData &data, SearchNode &search,
                                const LifeState &state, Placement p,
                                ContactType contactType) {
  const CatalystData &catalyst = data.catalysts[p.catalystIx];

  // Check that the type of catalyst contact matches the
  // active pattern neighbour count for an contact to occur
  if (catalyst.contactType != contactType)
    return PlacementValidity::INVALID_CONTACT;

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

  for (unsigned g = search.lookahead.gen; g < problem.gen; g++) {
    if constexpr (debug) std::cout << "Gen " << g << " state: " << current << std::endl;

    LifeState currentCount1(UNINITIALIZED), currentCount2(UNINITIALIZED),
        currentCountM(UNINITIALIZED);
    current.InteractionCounts(currentCount1, currentCount2,
                                          currentCountM);
    LifeState newContactPoints = (currentCount1 & ~currentHistory1) |
                                 (currentCount2 & ~currentHistory2) |
                                 (currentCountM & ~currentHistoryM);

    newContactPoints &= problem.LightCone(g);

    if constexpr (debug) std::cout << "Gen " << g << " contactPoints: " << newContactPoints << std::endl;

    LifeState remainingCells = newContactPoints;
    for (auto cell = remainingCells.FirstOn(); cell != std::make_pair(-1, -1);
         remainingCells.Erase(cell), cell = remainingCells.FirstOn()) {

      ContactType contactType =
          currentCount1.Get(cell)
              ? CONTACT1
              : (currentCount2.Get(cell) ? CONTACT2 : CONTACTM);

      for (unsigned i = 0; i < data.catalysts.size(); i++) {
        const CatalystData &catalyst = data.catalysts[i];

        Placement p = {cell, i};

        if (search.constraints[i].tried.Get(cell) || search.constraints[i].knownUnplaceable.Get(cell)) {
          continue;
        }

        if (catalyst.contactType == contactType && search.constraints[i].knownPlaceable.Get(cell)) {
          result.push_back(p);
          continue;
        }

        PlacementValidity validity =
            TestPlacement(data, search, current, p, contactType);

        switch (validity) {
        case PlacementValidity::VALID:
          search.constraints[i].knownPlaceable.Set(cell);
          result.push_back(p);
          break;
        case PlacementValidity::FAILED_CONTACT:
          search.constraints[i].knownUnplaceable.Set(cell);
          break;
        case PlacementValidity::FAILED_ELSEWHERE:
          break;
        case PlacementValidity::INVALID_CONTACT:
          // Might interact correctly in a future generation
          break;
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

std::vector<Perturbation>
CollatePerturbations(const SearchParams &params, const SearchData &data,
                     SearchNode &search, std::vector<Placement> &placements) {

  // TODO placeholder, actually dedup the perturbations
  std::vector<Perturbation> result;
  for (auto &p : placements) {
    result.push_back({0, {p}});
  }
  return result;
}

void MakePlacement(const SearchParams &params, const SearchData &data,
                   SearchNode &search, const Perturbation &p) {
  const Placement &placement = p.placements[0];
  const CatalystData &catalystdata = data.catalysts[placement.catalystIx];

  const LifeState catalyst = catalystdata.state.Moved(placement.pos);

  search.config.numCatalysts++;
  search.config.state |= catalyst;
  search.config.catalysts |= catalyst;
  search.config.required |= catalystdata.required.Moved(placement.pos);
  search.config.placements.push_back(placement);
  search.config.targets.push_back(LifeTarget(catalystdata.state.Moved(placement.pos), catalystdata.halo.Moved(placement.pos)));

  search.lookahead.state |= catalyst;
  search.lookahead.missingTime.push_back(0);

  search.history1 |= catalystdata.history1.Moved(placement.pos);
  search.history2 |= catalystdata.history2.Moved(placement.pos);
  search.historyM |= catalystdata.historyM.Moved(placement.pos);

  for (unsigned t = 0; t < data.catalysts.size(); t++) {
    search.constraints[t].tried |= data.collisionMasks[placement.catalystIx * data.catalysts.size() + t].Moved(placement.pos.first, placement.pos.second);
  }
}

void ResetLightcone(const SearchParams &params, const SearchData &data,
                    SearchNode &search, const Perturbation &p) {
  const Placement &placement = p.placements[0];

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
    tooClose = tooClose.ZOI(); // Is it faster to just recompute the `NZOIAround`?
  }

  // Now invalidate every placement that isn't safe
  // (Doesn't reset the `tried` field)
  for (unsigned i = 0; i < data.catalysts.size(); i++) {
    search.constraints[i].knownPlaceable &= safeContacts;
    search.constraints[i].knownUnplaceable &= safeContacts;
  }
}

void RunSearch(const SearchParams &params, const SearchData &data,
               SearchNode &search);

void BranchPerturbation(const SearchParams &params, const SearchData &data,
                        SearchNode &search, const Perturbation &p) {

  if constexpr (debug) std::cout << "Branching: " << p.placements[0].catalystIx << std::endl;

  for (auto &placement : p.placements) {
    search.constraints[placement.catalystIx].tried.Set(placement.pos);
  }

  SearchNode newSearch = search;

  MakePlacement(params, data, newSearch, p);

  if constexpr (debug) if (params.hasOracle && !(newSearch.config.state & ~params.oracle).IsEmpty())
    return;

  ResetLightcone(params, data, newSearch, p);

  // TODO?

  RunSearch(params, data, newSearch);
}

void RunSearch(const SearchParams &params, const SearchData &data,
               SearchNode &search) {

  if constexpr (debug) std::cout << "Starting node: " << search.lookahead.state << std::endl;

  TryAdvance(params, data, search);

  Problem problem = DetermineProblem(params, data, search.config, search.lookahead);

  if constexpr (print_progress) {
    static unsigned counter = 0;
    counter++;
    if (counter == 100000) {
      std::cout << "Current configuration: " << search.config.state << std::endl;
      LifeState problemGen = search.lookahead.state;
      problemGen.Step(problem.gen - search.lookahead.gen);
      std::cout << "Current problem: " << problem << std::endl << LifeHistoryState(problemGen, LifeState(), LifeState::Cell(problem.cell)) << std::endl;

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

  if (problem.type == ProblemType::NONE) {
    // TODO report properly!
    std::cout << "Winner: " << search.config.state << std::endl;
    if constexpr (debug) {
      std::cout << "Required: " << LifeHistoryState(search.config.state, LifeState(), search.config.required) << std::endl;
      for (auto &p : search.config.placements) {
        std::cout << "Placement: " << p.catalystIx << " at (" << p.pos.first << ", " << p.pos.second << ")" << std::endl;
      }
    }
    return;
  }

  if (search.config.numCatalysts == params.maxCatalysts)
    return;

  std::vector<Placement> placements =
      CollectPlacements(params, data, search, problem);

  std::vector<Perturbation> perturbations =
      CollatePerturbations(params, data, search, placements);

  for (auto &perturbation : perturbations) {
    BranchPerturbation(params, data, search, perturbation);
  }
}

std::vector<LifeState> CalculateCollisionMasks(const std::vector<CatalystData> &catalysts) {
  unsigned count = catalysts.size();
  std::vector<LifeState> result(count*count);
  for (unsigned s = 0; s < count; s++) {
      for (unsigned t = 0; t < count; t++) {
        result[s * count + t] = catalysts[s].CollisionMask(catalysts[t]);
      }
  }
  return result;
}

int main(int, char *argv[]) {
  auto toml = toml::parse(argv[1]);
  SearchParams params = SearchParams::FromToml(toml);

  std::vector<CatalystData> catalystdata;
  for (auto &c : params.catalysts) {
    auto newdata = CatalystData::FromParamsSymmetrically(c);
    catalystdata.insert(catalystdata.end(), newdata.begin(), newdata.end()); // Why is C++ like this
  }

  std::vector<LifeState> masks = CalculateCollisionMasks(catalystdata);

  SearchData data = {catalystdata, masks};

  SearchNode search(params, data);

  RunSearch(params, data, search);
}
