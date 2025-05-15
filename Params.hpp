#pragma once

#include "assert.h"

#include "toml.hpp"

#include "LifeAPI/LifeAPI.hpp"
#include "LifeAPI/LifeHistory.hpp"

enum class FilterType { EXACT, EVER, MATCH };

struct Filter {
  LifeState mask;
  LifeState state;
  std::pair<unsigned, unsigned> range;
  FilterType type;
};

// TODO
struct Forbidden {
  LifeState mask;
  LifeState state;
};

struct CatalystParams {
  LifeState state;
  LifeHistory required;
  std::vector<LifeHistory> approaches;
  std::vector<LifeHistory> forbiddens;
  std::vector<LifeState> soups;

  unsigned minRecoveryTime;
  unsigned maxRecoveryTime;

  bool transparent;
  bool limited;

  static CatalystParams FromToml(toml::value &toml);
};

CatalystParams CatalystParams::FromToml(toml::value &toml) {
  std::string rle = toml::find<std::string>(toml, "rle");
  LifeState state = LifeState::Parse(rle);

  bool transparent = toml::find_or(toml, "transparent", false);
  bool limited = toml::find_or(toml, "limited", false);

  LifeHistory required;
  std::vector<LifeHistory> approaches;

  if (!transparent) {
    // Not required if the catalyst is transparent
    std::string requiredrle = toml::find<std::string>(toml, "required");
    required = LifeHistory::Parse(requiredrle);

    if(toml.contains("approach")) {
      std::string rle = toml::find<std::string>(toml, "approach");
      approaches.push_back(LifeHistory::Parse(rle));
    }
    if(toml.contains("approaches")) {
      std::vector<std::string> approachrles = toml::find<std::vector<std::string>>(toml, "approaches");
      for (std::string &rle : approachrles)
        approaches.push_back(LifeHistory::Parse(rle));
    }
  }

  std::vector<LifeHistory> forbiddens;
  if (toml.contains("forbidden")) {
      std::vector<std::string> forbiddenrles = toml::find<std::vector<std::string>>(toml, "forbidden");
      for (std::string &rle : forbiddenrles)
        forbiddens.push_back(LifeHistory::Parse(rle));
  }

  std::vector<LifeState> soups;
  if(toml.contains("soup")) {
    std::string rle = toml::find<std::string>(toml, "soup");
    soups.push_back(LifeState::Parse(rle));
  }
  if(toml.contains("soups")) {
    std::vector<std::string> souprles = toml::find<std::vector<std::string>>(toml, "soups");
    for (std::string &rle : souprles)
      soups.push_back(LifeState::Parse(rle));
  }

  std::vector<int> recoveryRange =
      toml::find_or<std::vector<int>>(toml, "recovery-range", {0, 100});
  unsigned minRecoveryTime = recoveryRange[0];
  unsigned maxRecoveryTime = recoveryRange[1];

  return {state, required, approaches, forbiddens, soups,
          minRecoveryTime, maxRecoveryTime, transparent, limited};
}

struct SearchParams {
  std::vector<CatalystParams> catalysts;

  LifeHistory state;

  unsigned maxCatalysts;
  unsigned maxTransparent;
  unsigned maxLimited;

  unsigned minStableTime;

  // TODO are we actually checking all of these?
  unsigned minFirstActiveGen;
  unsigned maxFirstActiveGen;
  unsigned minActiveWindowGens;
  unsigned maxActiveWindowGens;
  unsigned maxStationaryTime;
  unsigned maxStationaryCount;

  bool continueAfterSuccess;
  unsigned maxStableTime;

  bool useBloomFilter;

  bool hasFilter;
  std::vector<Filter> filters;

  bool hasForbidden;
  std::vector<Forbidden> forbiddens;

  bool printSummary;

  std::string outputFile;

  bool debug;
  bool hasOracle;
  LifeState oracle;

  static SearchParams FromToml(toml::value &toml);
};

SearchParams SearchParams::FromToml(toml::value &toml) {
  SearchParams params;

  std::string rle = toml::find<std::string>(toml, "pattern");
  params.state = LifeHistory::Parse(rle);

  params.maxCatalysts = toml::find_or(toml, "max-catalysts", 100);
  params.maxTransparent = toml::find_or(toml, "max-transparent", 0);
  params.maxLimited = toml::find_or(toml, "max-limited", 1);
  params.minStableTime = toml::find_or(toml, "min-stable-time", 8);
  params.continueAfterSuccess = toml::find_or(toml, "continue-after-success", true);
  params.maxStableTime = toml::find_or(toml, "max-stable-time", std::max(params.minStableTime, (unsigned)64));

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

  params.printSummary = toml::find_or(toml, "print-summary", false);
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
      LifeHistory pat = LifeHistory::Parse(rle);
      pat.Move(filterCenterVec[0], filterCenterVec[1]);

      int filterGen = toml::find_or(f, "filter-gen", -1);
      std::pair<int, int> filterRange;
      if (filterGen != -1) {
        filterRange = {filterGen, filterGen};
      } else {
        std::vector<int> filterRangeVec = toml::find_or<std::vector<int>>(f, "filter-range", {-1, -1});
        filterRange = {filterRangeVec[0], filterRangeVec[1]};
      }

      FilterType filterType;
      std::string filterTypeStr =
          toml::find_or<std::string>(f, "filter-type", "EVER");
      if (filterTypeStr == "EXACT") {
        filterType = FilterType::EXACT;
      } else if (filterTypeStr == "EVER") {
        filterType = FilterType::EVER;
      } else if (filterTypeStr == "MATCH") {
        filterType = FilterType::MATCH;
      }

      params.filters.push_back({pat.marked, pat.state, filterRange, filterType});
    }
  } else {
    params.hasFilter = false;
  }

  params.debug = toml::find_or(toml, "debug", false);

  if (toml.contains("oracle")) {
    std::string oraclerle = toml::find<std::string>(toml, "oracle");
    LifeHistory oracle = LifeHistory::Parse(oraclerle);

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
