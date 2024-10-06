#pragma once

#include "assert.h"

#include "toml.hpp"

#include "LifeAPI/LifeAPI.hpp"
#include "LifeAPI/LifeHistoryState.hpp"

enum class FilterType { EXACT, EVER };

// TODO
struct Filter {
  LifeState mask;
  LifeState state;
  unsigned gen;
  FilterType type;
};

// TODO
struct Forbidden {
  LifeState mask;
  LifeState state;
};

struct CatalystParams {
  LifeState state;
  LifeHistoryState required;
  std::vector<LifeHistoryState> approaches;
  std::vector<LifeHistoryState> forbiddens;
  std::vector<LifeState> soups;

  unsigned minRecoveryTime;
  unsigned maxRecoveryTime;

  bool transparent;

  static CatalystParams FromToml(toml::value &toml);
};

CatalystParams CatalystParams::FromToml(toml::value &toml) {
  std::string rle = toml::find<std::string>(toml, "rle");
  LifeState state = LifeState::Parse(rle);

  bool transparent = toml::find_or(toml, "transparent", false);

  LifeHistoryState required;
  std::vector<LifeHistoryState> approaches;

  if (!transparent) {
    // Not required if the catalyst is transparent
    std::string requiredrle = toml::find<std::string>(toml, "required");
    required = LifeHistoryState::Parse(requiredrle);

    if(toml.contains("approach")) {
      std::string rle = toml::find<std::string>(toml, "approach");
      approaches.push_back(LifeHistoryState::Parse(rle));
    }
    if(toml.contains("approaches")) {
      std::vector<std::string> approachrles = toml::find<std::vector<std::string>>(toml, "approaches");
      for (std::string &rle : approachrles)
        approaches.push_back(LifeHistoryState::Parse(rle));
    }
  }

  std::vector<LifeHistoryState> forbiddens;
  if (toml.contains("forbidden")) {
      std::vector<std::string> forbiddenrles = toml::find<std::vector<std::string>>(toml, "forbidden");
      for (std::string &rle : forbiddenrles)
        forbiddens.push_back(LifeHistoryState::Parse(rle));
  }

  std::vector<int> recoveryRange =
      toml::find_or<std::vector<int>>(toml, "recovery-range", {0, 100});
  unsigned minRecoveryTime = recoveryRange[0];
  unsigned maxRecoveryTime = recoveryRange[1];

  return {state, required, approaches, forbiddens, std::vector<LifeState>(),
      minRecoveryTime, maxRecoveryTime, transparent};
}

struct SearchParams {
  std::vector<CatalystParams> catalysts;

  LifeHistoryState state;

  unsigned maxCatalysts;
  unsigned maxTransparent;

  unsigned minStableTime;

  // TODO are we actually checking all of these?
  unsigned minFirstActiveGen;
  unsigned maxFirstActiveGen;
  unsigned minActiveWindowGens;
  unsigned maxActiveWindowGens;
  unsigned maxStationaryTime;
  unsigned maxStationaryCount;

  bool useBloomFilter;

  bool hasFilter;
  std::vector<Filter> filters;

  bool hasForbidden;
  std::vector<Forbidden> forbiddens;

  std::string outputFile;

  bool debug;
  bool hasOracle;
  LifeState oracle;

  static SearchParams FromToml(toml::value &toml);
};

SearchParams SearchParams::FromToml(toml::value &toml) {
  SearchParams params;

  std::string rle = toml::find<std::string>(toml, "pattern");
  params.state = LifeHistoryState::Parse(rle);

  params.maxCatalysts = toml::find_or(toml, "max-catalysts", 100);
  params.maxTransparent = toml::find_or(toml, "max-transparent", 0);
  params.minStableTime = toml::find_or(toml, "min-stable-time", 8);

  std::vector<int> firstRange =
      toml::find_or<std::vector<int>>(toml, "first-active-range", {0, 100});
  params.minFirstActiveGen = firstRange[0];
  params.maxFirstActiveGen = firstRange[1];

  std::vector<int> windowRange = toml::find_or<std::vector<int>>(toml, "active-window-range", {0, 500});
  params.minActiveWindowGens = windowRange[0];
  params.maxActiveWindowGens = windowRange[1];

  params.maxStationaryTime = toml::find_or(toml, "max-stationary-time", 0);
  params.maxStationaryCount = toml::find_or(toml, "max-stationary-count", 0);

  params.useBloomFilter = toml::find_or(toml, "use-bloom-filter", false);

  params.outputFile = toml::find_or(toml, "output-file", "lightcone-output.rle");

  std::vector<int> patternCenterVec = toml::find_or<std::vector<int>>(toml, "pattern-center", {0, 0});
  std::pair<int, int> patternCenter = {-patternCenterVec[0], -patternCenterVec[1]};
  params.state.Move(patternCenter);

  if (toml.contains("catalyst")) {
    auto catalysts = toml::find<std::vector<toml::value>>(toml, "catalyst");
    for (auto &c : catalysts) {
      params.catalysts.push_back(CatalystParams::FromToml(c));
    }
  }

  if (toml.contains("filter")) {
    params.hasFilter = true;

    auto filters = toml::find<std::vector<toml::value>>(toml, "filter");
    for (auto &f : filters) {
      std::string rle = toml::find_or<std::string>(f, "filter", "");
      std::vector<int> filterCenterVec =
          toml::find_or<std::vector<int>>(f, "filter-pos", {0, 0});
      LifeHistoryState pat = LifeHistoryState::Parse(rle);
      unsigned filterGen = toml::find_or(f, "filter-gen", -1);

      FilterType filterType;
      std::string filterTypeStr =
          toml::find_or<std::string>(f, "filter-type", "EVER");
      if (filterTypeStr == "EXACT") {
        filterType = FilterType::EXACT;
      }
      if (filterTypeStr == "EVER") {
        filterType = FilterType::EVER;
      }

      pat.Move(filterCenterVec[0], filterCenterVec[1]);

      params.filters.push_back({pat.marked, pat.state, filterGen, filterType});
    }
  } else {
    params.hasFilter = false;
  }

  params.debug = toml::find_or(toml, "debug", false);

  if (toml.contains("oracle")) {
    std::string oraclerle = toml::find<std::string>(toml, "oracle");
    LifeHistoryState oracle = LifeHistoryState::Parse(oraclerle);

    std::vector<int> oracleCenterVec = toml::find_or<std::vector<int>>(toml, "oracle-center", {0, 0});
    std::pair<int, int> oracleCenter = {-oracleCenterVec[0], -oracleCenterVec[1]};
    oracle.Move(oracleCenter);

    params.hasOracle = true;
    params.oracle = oracle.state;
  } else {
    params.hasOracle = false;
  }

  return params;
}
