#include "CbcSolverStatistics.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

std::string toLower(std::string value)
{
  std::transform(value.begin(), value.end(), value.begin(),
    [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  return value;
}

std::string stripExtension(const std::string &filename)
{
  std::string base = filename;
  if (base.size() > 3 && base.compare(base.size() - 3, 3, ".gz") == 0)
    base = base.substr(0, base.size() - 3);
  auto endsWith = [&](const std::string &ext) {
    return base.size() >= ext.size() &&
      base.compare(base.size() - ext.size(), ext.size(), ext) == 0;
  };
  if (endsWith(".mps"))
    return base.substr(0, base.size() - 4);
  if (endsWith(".lp"))
    return base.substr(0, base.size() - 3);
  return filename;
}

std::string stripPath(const std::string &value)
{
  std::string::size_type pos = value.find_last_of("/\\");
  if (pos == std::string::npos)
    return value;
  return value.substr(pos + 1);
}

std::string buildRuntimeOptions(const std::deque<std::string> &tokens)
{
  std::ostringstream stream;
  bool first = true;
  for (const std::string &token : tokens) {
    if (token.empty())
      continue;
    if (token == "cbc" || token == "clp")
      continue;
    std::string lower = toLower(token);
    if (lower.find(".mps") != std::string::npos ||
      lower.find(".gz") != std::string::npos)
      continue;
    if (lower.rfind("-writestat", 0) == 0)
      break;
    if (!first)
      stream << ' ';
    stream << token;
    first = false;
  }
  return stream.str();
}

/**
 * Convert a generator/heuristic name to a safe column-name fragment.
 * Lowercases, replaces runs of non-alphanumeric characters with a single
 * underscore, and strips leading/trailing underscores.
 */
std::string sanitizeName(const std::string &name)
{
  std::string result;
  result.reserve(name.size());
  bool lastWasUnderscore = true; // suppress leading underscore
  for (unsigned char ch : name) {
    if (std::isalnum(ch)) {
      result += static_cast<char>(std::tolower(ch));
      lastWasUnderscore = false;
    } else {
      if (!lastWasUnderscore)
        result += '_';
      lastWasUnderscore = true;
    }
  }
  // strip trailing underscore
  if (!result.empty() && result.back() == '_')
    result.pop_back();
  return result;
}

std::string formatDouble(double value, int precision,
  std::ios_base::fmtflags floatField = std::ios_base::fmtflags(0))
{
  std::ostringstream out;
  if (floatField != std::ios_base::fmtflags(0))
    out.setf(floatField, std::ios_base::floatfield);
  std::streamsize oldPrecision = out.precision();
  out << std::setprecision(precision) << value;
  out.precision(oldPrecision);
  return out.str();
}

} // namespace

// ---------------------------------------------------------------------------
// Canonical lists
// ---------------------------------------------------------------------------

const std::vector<std::string> &CbcSolverStatistics::knownCutGenerators()
{
  // Names exactly as passed to CbcModel::addCutGenerator() in
  // CbcSolverCutSetup.cpp.  Add new names here when new generators are
  // introduced so that the CSV column set stays stable.
  static const std::vector<std::string> kList = {
    "Probing",
    "Gomory",
    "GomoryL1",
    "GomoryL2",
    "Gomory(2)",
    "Knapsack",
    "Reduce-and-split",
    "Reduce-and-split(2)",
    "Clique",
    "OddWheel",
    "MixedIntegerRounding2",
    "FlowCover",
    "TwoMirCuts",
    "TwoMirCutsL1",
    "TwoMirCutsL2",
    "LiftAndProject",
    "ResidualCapacity",
    "ZeroHalf",
    "Stored",
  };
  return kList;
}

const std::vector<std::string> &CbcSolverStatistics::knownHeuristics()
{
  // Names exactly as passed to setHeuristicName() in CbcSolverHeuristics.cpp
  // and CbcSolver.cpp.  Names that sanitize to the same string are treated
  // as the same heuristic (e.g. "feasibility pump" and "Feasibility pump").
  static const std::vector<std::string> kList = {
    "feasibility pump",
    "rounding",
    "combine solutions",
    "greedy cover",
    "greedy equality",
    "random rounding",
    "dynamic pass thru",
    "linked",
    "Partial solution given",
    "FeasibilityJump",
    "Dantzig-Wolfe-expansion",
    "RINS",
    "RENS",
    "RENSdj",
    "RENSub",
    "VND",
    "Naive",
    "DiveAny",
    "DiveCoefficient",
    "DiveFractional",
    "DiveGuided",
    "DiveLineSearch",
    "DivePseudoCost",
    "DiveVectorLength",
    "Multiple root solvers",
  };
  return kList;
}

// ---------------------------------------------------------------------------
// writeCsv
// ---------------------------------------------------------------------------

bool CbcSolverStatistics::writeCsv(CbcParameters &parameters,
  const std::string &outFileName,
  const std::deque<std::string> &inputQueue) const
{
  if (outFileName.empty())
    return false;

  // ------------------------------------------------------------------
  // 1. Build cut and heuristic column lists from the canonical sets only.
  //    This guarantees the header is identical across all runs regardless
  //    of which generators/heuristics were active, so we can always append
  //    without ever needing to rewrite existing rows.
  // ------------------------------------------------------------------
  std::vector<std::string> cutSanitized;
  for (const auto &n : knownCutGenerators())
    cutSanitized.push_back(sanitizeName(n));

  std::vector<std::string> heurSanitized;
  for (const auto &n : knownHeuristics())
    heurSanitized.push_back(sanitizeName(n));

  // ------------------------------------------------------------------
  // 2. Build maps: sanitized-name -> aggregated stats (accumulate
  //    entries with the same sanitized name, e.g. "feasibility pump"
  //    and "Feasibility pump").
  // ------------------------------------------------------------------
  std::unordered_map<std::string, CutGeneratorStats> cutMap;
  for (const auto &cs : cutStats) {
    std::string san = sanitizeName(cs.name);
    auto &acc = cutMap[san];
    if (acc.name.empty())
      acc.name = cs.name; // first entry wins
    acc.nCuts += cs.nCuts;
    acc.nCalls += cs.nCalls;
    acc.time += cs.time;
    acc.nColumnCuts += cs.nColumnCuts;
    acc.nElements += cs.nElements;
    if (cs.minNz >= 0 && (acc.minNz < 0 || cs.minNz < acc.minNz))
      acc.minNz = cs.minNz;
    if (cs.maxNz > acc.maxNz)
      acc.maxNz = cs.maxNz;
    if (cs.minDepth >= 0 && (acc.minDepth < 0 || cs.minDepth < acc.minDepth))
      acc.minDepth = cs.minDepth;
    if (cs.maxDepth > acc.maxDepth)
      acc.maxDepth = cs.maxDepth;
  }

  std::unordered_map<std::string, HeuristicStats> heurMap;
  for (const auto &hs : heuristicStats) {
    std::string san = sanitizeName(hs.name);
    auto &acc = heurMap[san];
    if (acc.name.empty())
      acc.name = hs.name; // first entry wins
    acc.nExecutions += hs.nExecutions;
    acc.totalTime += hs.totalTime;
    acc.nSolutions += hs.nSolutions;
    if (hs.minDepth >= 0 && (acc.minDepth < 0 || hs.minDepth < acc.minDepth))
      acc.minDepth = hs.minDepth;
    if (hs.maxDepth > acc.maxDepth)
      acc.maxDepth = hs.maxDepth;
  }

  // ------------------------------------------------------------------
  // 3. Build the header string.
  //    For each cut:  cut_<san>_cuts, cut_<san>_calls, cut_<san>_time,
  //                   cut_<san>_elements, cut_<san>_minnz, cut_<san>_maxnz,
  //                   cut_<san>_avgnz, cut_<san>_minDepth, cut_<san>_maxDepth
  //    For each heur: heur_<san>_execs, heur_<san>_time, heur_<san>_sols,
  //                   heur_<san>_minDepth, heur_<san>_maxDepth
  // ------------------------------------------------------------------
  std::ostringstream headerStream;
  headerStream << "Name,result,integer_feasible,time,sys,elapsed,objective,continuous,"
               << "lp_seconds,tightened,cut_time,"
               << "nodes,iterations,rows,columns,processed_rows,"
               << "processed_columns,cgraph_time,cgraph_density";
  for (const auto &san : cutSanitized) {
    headerStream << ",cut_" << san << "_cuts"
                 << ",cut_" << san << "_calls"
                 << ",cut_" << san << "_time"
                 << ",cut_" << san << "_elements"
                 << ",cut_" << san << "_minnz"
                 << ",cut_" << san << "_maxnz"
                 << ",cut_" << san << "_avgnz"
                 << ",cut_" << san << "_minDepth"
                 << ",cut_" << san << "_maxDepth";
  }
  for (const auto &san : heurSanitized) {
    headerStream << ",heur_" << san << "_execs"
                 << ",heur_" << san << "_time"
                 << ",heur_" << san << "_sols"
                 << ",heur_" << san << "_minDepth"
                 << ",heur_" << san << "_maxDepth";
  }
  headerStream << ",runtime_options";
  const std::string headerLine = headerStream.str();

  // ------------------------------------------------------------------
  // 4. Build the new data line.
  // ------------------------------------------------------------------
  std::string inputFileName = parameters[CbcParam::IMPORTFILE]->strVal();
  const std::string problemName = stripExtension(stripPath(inputFileName));
  const std::string runtimeOptions = buildRuntimeOptions(inputQueue);

  std::ostringstream dataStream;
  dataStream << problemName << ',' << result << ','
             << (integer_feasible ? "1" : "0") << ','
             << formatDouble(seconds, 2, std::ios_base::fixed) << ','
             << formatDouble(sys_seconds, 2, std::ios_base::fixed) << ','
             << formatDouble(elapsed_seconds, 2, std::ios_base::fixed) << ','
             << formatDouble(obj, 16) << ','
             << formatDouble(continuous, 6) << ','
             << formatDouble(lp_seconds, 2, std::ios_base::fixed) << ','
             << formatDouble(tighter, 6) << ','
             << formatDouble(cut_time, 2, std::ios_base::fixed) << ','
             << nodes << ',' << iterations << ',' << nrows << ',' << ncols
             << ',' << nprocessedrows << ',' << nprocessedcols
             << ',' << formatDouble(cgraph_time, 2, std::ios_base::fixed)
             << ',' << formatDouble(cgraph_density, 6);

  for (const auto &san : cutSanitized) {
    const auto it = cutMap.find(san);
    if (it != cutMap.end()) {
      const auto &cs = it->second;
      double avgNz = (cs.nCuts > 0) ? static_cast<double>(cs.nElements) / cs.nCuts : 0.0;
      dataStream << ',' << cs.nCuts
                 << ',' << cs.nCalls
                 << ',' << formatDouble(cs.time, 4, std::ios_base::fixed)
                 << ',' << cs.nElements
                 << ',' << cs.minNz
                 << ',' << cs.maxNz
                 << ',' << formatDouble(avgNz, 2, std::ios_base::fixed)
                 << ',' << cs.minDepth
                 << ',' << cs.maxDepth;
    } else {
      dataStream << ",0,0,0.0000,0,-1,-1,0.00,-1,-1";
    }
  }
  for (const auto &san : heurSanitized) {
    const auto it = heurMap.find(san);
    if (it != heurMap.end()) {
      const auto &hs = it->second;
      dataStream << ',' << hs.nExecutions
                 << ',' << formatDouble(hs.totalTime, 4, std::ios_base::fixed)
                 << ',' << hs.nSolutions
                 << ',' << hs.minDepth
                 << ',' << hs.maxDepth;
    } else {
      dataStream << ",0,0.0000,0,-1,-1";
    }
  }
  dataStream << ',' << runtimeOptions;
  const std::string newDataLine = dataStream.str();

  // ------------------------------------------------------------------
  // 5. Write to file: write header only when the file is new or empty,
  //    then always append the data line.
  // ------------------------------------------------------------------
  bool hasHeader = false;
  {
    std::ifstream in(outFileName.c_str());
    std::string line;
    if (std::getline(in, line) && !line.empty())
      hasHeader = true;
  }

  std::ofstream file(outFileName.c_str(), std::ios::out | std::ios::app);
  if (!file.is_open())
    return false;
  if (!hasHeader)
    file << headerLine << '\n';
  file << newDataLine << std::endl;
  return true;
}
