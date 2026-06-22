#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include <map>
#include <set>
#include "CglLKCI.hpp"
#include "OsiRowCut.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinKnapsackRow.hpp"

CglLKCI::CglLKCI()
    : CglCutGenerator(),
      beam_width_(2),
      random_starts_(1),
      local_rounds_(2),
      move_pool_(8),
      seed_(0),
      tol_(1e-7),
      max_cover_rows_(-1),
      max_row_nnz_(-1) {}

CglLKCI::CglLKCI(const CglLKCI& rhs)
    : CglCutGenerator(rhs),
      beam_width_(rhs.beam_width_),
      random_starts_(rhs.random_starts_),
      local_rounds_(rhs.local_rounds_),
      move_pool_(rhs.move_pool_),
      seed_(rhs.seed_),
      tol_(rhs.tol_),
      max_cover_rows_(rhs.max_cover_rows_),
      max_row_nnz_(rhs.max_row_nnz_) {}

CglLKCI::~CglLKCI() {}

CglLKCI& CglLKCI::operator=(const CglLKCI& rhs) {
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    beam_width_ = rhs.beam_width_;
    random_starts_ = rhs.random_starts_;
    local_rounds_ = rhs.local_rounds_;
    move_pool_ = rhs.move_pool_;
    seed_ = rhs.seed_;
    tol_ = rhs.tol_;
    max_cover_rows_ = rhs.max_cover_rows_;
    max_row_nnz_ = rhs.max_row_nnz_;
  }
  return *this;
}

CglCutGenerator* CglLKCI::clone() const {
  return new CglLKCI(*this);
}

std::string CglLKCI::generateCpp(FILE* fp) {
  return "";
}

double CglLKCI::h_value(long long r, long long D, const std::vector<long long>& T, const std::vector<long long>& prefix_min) {
  int m = T.size() - 1;
  if (m == 0) {
    return (double)r;
  }
  long long z = r + D;
  auto it = std::lower_bound(T.begin(), T.end(), z);
  int p = std::distance(T.begin(), it);

  double term1 = (p <= m) ? (double)(p - 1) * D : std::numeric_limits<double>::infinity();
  int q = std::min(p - 1, m);
  if (q < 0) q = 0;
  double term2 = (double)r + (double)prefix_min[q];
  return std::min(term1, term2);
}

std::vector<double> CglLKCI::h_values(const std::vector<long long>& r, long long D, const std::vector<long long>& b_desc) {
  int m = b_desc.size();
  int n = r.size();
  std::vector<double> out(n);
  if (m == 0) {
    for (int i = 0; i < n; ++i) {
      out[i] = (double)r[i];
    }
    return out;
  }
  std::vector<long long> T(m + 1);
  T[0] = 0;
  for (int i = 0; i < m; ++i) {
    T[i + 1] = T[i] + b_desc[i];
  }
  std::vector<long long> prefix_min(m + 1);
  prefix_min[0] = 0;
  for (int i = 1; i <= m; ++i) {
    long long val = i * D - T[i];
    prefix_min[i] = std::min(prefix_min[i - 1], val);
  }
  for (int i = 0; i < n; ++i) {
    out[i] = h_value(r[i], D, T, prefix_min);
  }
  return out;
}

bool CglLKCI::evaluate_lkci(
    const std::vector<long long>& a,
    long long d,
    const std::vector<double>& x,
    const std::vector<bool>& S_mask,
    LKCICandidate& cand,
    double tol,
    bool require_positive)
{
  int n = a.size();
  long long S_sum = 0;
  std::vector<int> S_idx;
  for (int i = 0; i < n; ++i) {
    if (S_mask[i]) {
      S_sum += a[i];
      S_idx.push_back(i);
    }
  }
  long long D = d - S_sum;
  if (D <= 0) return false;

  std::vector<int> R_idx;
  long long R_sum = 0;
  long long R_min = std::numeric_limits<long long>::max();
  for (int i = 0; i < n; ++i) {
    if (!S_mask[i] && a[i] <= D) {
      R_idx.push_back(i);
      R_sum += a[i];
      if (a[i] < R_min) {
        R_min = a[i];
      }
    }
  }

  if (R_idx.empty()) return false;
  if (R_sum < D) return false;
  if (R_sum - R_min < D) return false;

  std::vector<int> L_idx;
  std::vector<long long> b_desc;
  for (int i = 0; i < n; ++i) {
    if (!S_mask[i] && a[i] > D) {
      L_idx.push_back(i);
      b_desc.push_back(a[i]);
    }
  }
  std::sort(b_desc.begin(), b_desc.end(), std::greater<long long>());

  double lhs = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!S_mask[i]) {
      double coeff = std::min((double)a[i], (double)D);
      lhs += coeff * x[i];
    }
  }

  std::vector<double> gamma(n, 0.0);
  double rhs = (double)D;
  if (!S_idx.empty()) {
    std::vector<long long> r_S(S_idx.size());
    for (size_t i = 0; i < S_idx.size(); ++i) {
      r_S[i] = a[S_idx[i]];
    }
    std::vector<double> gamma_S = h_values(r_S, D, b_desc);
    for (size_t i = 0; i < S_idx.size(); ++i) {
      int idx = S_idx[i];
      gamma[idx] = gamma_S[i];
      rhs += gamma_S[i] * (1.0 - x[idx]);
    }
  }

  double delta = rhs - lhs;
  if (require_positive && delta <= tol) {
    return false;
  }

  cand.S = S_idx;
  cand.D = D;
  cand.delta = delta;
  cand.L = L_idx;
  cand.R = R_idx;
  cand.gamma = gamma;
  cand.lhs = lhs;
  cand.rhs = rhs;
  return true;
}

std::vector<int> CglLKCI::topk_indices(const std::vector<double>& scores, const std::vector<bool>& mask, int k) {
  std::vector<int> idx;
  for (size_t i = 0; i < mask.size(); ++i) {
    if (mask[i]) {
      idx.push_back(i);
    }
  }
  if (idx.empty()) return idx;
  int n = idx.size();
  int num_keep = std::min(k, n);

  std::sort(idx.begin(), idx.end(), [&](int i, int j) {
    return scores[i] > scores[j];
  });

  if (idx.size() > (size_t)num_keep) {
    idx.resize(num_keep);
  }
  return idx;
}

void CglLKCI::add_prefix_candidates(
    const std::vector<long long>& a,
    long long d,
    const std::vector<double>& x,
    const std::vector<int>& order,
    const std::vector<long long>& targets,
    std::vector<BeamEntry>& beam,
    int beam_width,
    double tol)
{
  int n = a.size();
  std::vector<bool> mask(n, false);
  long long total = 0;

  std::vector<long long> filtered_targets;
  for (long long t : targets) {
    if (t > 0 && t < d) {
      filtered_targets.push_back(t);
    }
  }
  std::sort(filtered_targets.begin(), filtered_targets.end(), std::greater<long long>());

  size_t target_pos = 0;
  long long next_target = filtered_targets.empty() ? -1 : filtered_targets[0];

  auto consider = [&]() {
    LKCICandidate cand;
    if (evaluate_lkci(a, d, x, mask, cand, tol, false)) {
      BeamEntry entry{ cand.delta, mask, cand };
      beam.push_back(entry);
      std::sort(beam.begin(), beam.end());
      if (beam.size() > (size_t)beam_width) {
        beam.resize(beam_width);
      }
    }
  };

  consider();

  for (size_t t = 0; t < order.size(); ++t) {
    int j = order[t];
    long long aj = a[j];

    if (total + aj >= d) {
      continue;
    }

    mask[j] = true;
    total += aj;
    long long D = d - total;

    bool should_eval = (t < 32) || (((t + 1) & t) == 0);

    while (next_target != -1 && D <= next_target) {
      should_eval = true;
      target_pos++;
      next_target = (target_pos < filtered_targets.size()) ? filtered_targets[target_pos] : -1;
    }

    if (should_eval) {
      consider();
    }
  }
}

bool CglLKCI::run_exact_best_move(
    const std::vector<long long>& a,
    long long d,
    const std::vector<double>& x,
    const std::vector<long long>& targets_vec,
    std::vector<bool>& mask,
    LKCICandidate& cand,
    double tol,
    int move_pool,
    double& best_gain)
{
  int n = a.size();
  double Lx = 0.0;
  for (int j : cand.L) {
    Lx += x[j];
  }

  std::vector<long long> b_desc;
  for (int j : cand.L) {
    b_desc.push_back(a[j]);
  }
  std::sort(b_desc.begin(), b_desc.end(), std::greater<long long>());

  std::vector<double> h_all = h_values(a, cand.D, b_desc);

  std::vector<double> add_proxy(n);
  std::vector<bool> add_mask(n);
  for (int i = 0; i < n; ++i) {
    add_proxy[i] = std::min((double)a[i], (double)cand.D) * x[i] + h_all[i] * (1.0 - x[i]) + (double)a[i] * (Lx - 1.0);
    add_mask[i] = !mask[i] && (a[i] < cand.D);
  }

  std::vector<int> add_idx = topk_indices(add_proxy, add_mask, move_pool);

  std::vector<double> rem_proxy(n, -std::numeric_limits<double>::infinity());
  int num_in_mask = 0;
  for (int i = 0; i < n; ++i) {
    if (mask[i]) {
      rem_proxy[i] = -(h_all[i] * (1.0 - x[i]) + std::min((double)a[i], (double)cand.D) * x[i]);
      num_in_mask++;
    }
  }

  std::vector<int> rem_idx = topk_indices(rem_proxy, mask, std::min(move_pool, num_in_mask));

  LKCICandidate best_cand = cand;
  std::vector<bool> best_move_mask = mask;
  best_gain = 0.0;

  for (int j : add_idx) {
    std::vector<bool> m = mask;
    m[j] = true;
    LKCICandidate c;
    if (evaluate_lkci(a, d, x, m, c, tol, false)) {
      if (c.delta - cand.delta > best_gain + tol) {
        best_cand = c;
        best_move_mask = m;
        best_gain = c.delta - cand.delta;
      }
    }
  }

  for (int j : rem_idx) {
    std::vector<bool> m = mask;
    m[j] = false;
    LKCICandidate c;
    if (evaluate_lkci(a, d, x, m, c, tol, false)) {
      if (c.delta - cand.delta > best_gain + tol) {
        best_cand = c;
        best_move_mask = m;
        best_gain = c.delta - cand.delta;
      }
    }
  }

  if (best_gain <= tol && !add_idx.empty() && !rem_idx.empty()) {
    int max_r = std::min(8, (int)rem_idx.size());
    int max_j = std::min(8, (int)add_idx.size());
    long long current_sum_a = 0;
    for (int i = 0; i < n; ++i) {
      if (mask[i]) current_sum_a += a[i];
    }

    for (int ri = 0; ri < max_r; ++ri) {
      int r = rem_idx[ri];
      for (int ji = 0; ji < max_j; ++ji) {
        int j = add_idx[ji];
        if (j == r) continue;

        long long new_D = d - current_sum_a + a[r] - a[j];
        if (new_D <= 0) continue;

        std::vector<bool> m = mask;
        m[r] = false;
        m[j] = true;

        LKCICandidate c;
        if (evaluate_lkci(a, d, x, m, c, tol, false)) {
          if (c.delta - cand.delta > best_gain + tol) {
            best_cand = c;
            best_move_mask = m;
            best_gain = c.delta - cand.delta;
          }
        }
      }
    }
  }

  if (best_gain > tol) {
    mask = best_move_mask;
    cand = best_cand;
    return true;
  }
  return false;
}

bool CglLKCI::separate_lkci(
    const std::vector<long long>& a,
    long long d,
    const std::vector<double>& x,
    LKCICandidate& cand_out,
    int seed_offset)
{
  int n = a.size();
  std::mt19937 generator(seed_ + seed_offset);
  std::uniform_real_distribution<double> distribution(1e-10, 1.0 - 1e-10);

  std::set<long long> targets;
  std::vector<long long> pos;
  for (long long val : a) {
    if (val > 0 && val < d) {
      pos.push_back(val);
    }
  }
  if (!pos.empty()) {
    std::sort(pos.begin(), pos.end());
    int M = pos.size();
    for (int i = 1; i <= 19; ++i) {
      double f = 0.05 * i;
      int idx = (int)std::round(f * (M - 1));
      if (idx < 0) idx = 0;
      if (idx >= M) idx = M - 1;
      long long q = pos[idx];
      long long t = std::max(1LL, std::min(d - 1, q));
      targets.insert(t);
    }
  }
  std::vector<double> fixed_fracs = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05};
  for (double f : fixed_fracs) {
    long long t = std::max(1LL, std::min(d - 1, (long long)std::round((double)d * f)));
    targets.insert(t);
  }

  std::vector<long long> targets_vec(targets.begin(), targets.end());
  std::vector<std::vector<int>> orders;

  auto get_order = [&](const std::vector<double>& scores) {
    std::vector<int> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int i, int j) {
      return scores[i] > scores[j];
    });
    return ord;
  };

  // 1. a * (1.0 - x)
  {
    std::vector<double> s(n);
    for (int i = 0; i < n; ++i) s[i] = (double)a[i] * (1.0 - x[i]);
    orders.push_back(get_order(s));
  }
  // 2. a * x
  {
    std::vector<double> s(n);
    for (int i = 0; i < n; ++i) s[i] = (double)a[i] * x[i];
    orders.push_back(get_order(s));
  }
  // 3. a / (x + 0.05)
  {
    std::vector<double> s(n);
    for (int i = 0; i < n; ++i) s[i] = (double)a[i] / (x[i] + 0.05);
    orders.push_back(get_order(s));
  }
  // 4. a * abs(x - 0.5)
  {
    std::vector<double> s(n);
    for (int i = 0; i < n; ++i) s[i] = (double)a[i] * std::abs(x[i] - 0.5);
    orders.push_back(get_order(s));
  }
  // 5. a as float
  {
    std::vector<double> s(n);
    for (int i = 0; i < n; ++i) s[i] = (double)a[i];
    orders.push_back(get_order(s));
  }

  // 6. Gumbel random starts
  std::vector<double> base(n);
  for (int i = 0; i < n; ++i) {
    base[i] = (double)a[i] * (0.25 + std::max(x[i], 1.0 - x[i]));
  }
  for (int r = 0; r < random_starts_; ++r) {
    std::vector<double> s(n);
    for (int i = 0; i < n; ++i) {
      double u = distribution(generator);
      double gumbel = -std::log(-std::log(u));
      s[i] = std::log(base[i] + 1e-9) + 0.4 * gumbel;
    }
    orders.push_back(get_order(s));
  }

  std::vector<BeamEntry> beam;
  for (const auto& order : orders) {
    add_prefix_candidates(a, d, x, order, targets_vec, beam, beam_width_, tol_);
  }

  if (beam.empty()) {
    std::vector<double> scores(n);
    std::vector<bool> mask_a(n);
    for (int i = 0; i < n; ++i) {
      scores[i] = (double)a[i];
      mask_a[i] = (a[i] < d);
    }
    std::vector<int> cand_items = topk_indices(scores, mask_a, std::min(64, n));
    for (int j : cand_items) {
      std::vector<bool> mask(n, false);
      mask[j] = true;
      LKCICandidate cand;
      if (evaluate_lkci(a, d, x, mask, cand, tol_, false)) {
        beam.push_back(BeamEntry{ cand.delta, mask, cand });
      }
    }
    std::sort(beam.begin(), beam.end());
    if (beam.size() > (size_t)beam_width_) {
      beam.resize(beam_width_);
    }
  }

  if (beam.empty()) {
    return false;
  }

  LKCICandidate best = beam[0].cand;
  std::vector<bool> best_mask = beam[0].mask;

  for (const auto& entry : beam) {
    std::vector<bool> mask = entry.mask;
    LKCICandidate cand = entry.cand;

    for (int round = 0; round < local_rounds_; ++round) {
      double gain = 0.0;
      if (!run_exact_best_move(a, d, x, targets_vec, mask, cand, tol_, move_pool_, gain)) {
        break;
      }
    }

    if (cand.delta > best.delta) {
      best = cand;
      best_mask = mask;
    }
  }

  LKCICandidate final_cand;
  if (evaluate_lkci(a, d, x, best_mask, final_cand, tol_, true)) {
    cand_out = final_cand;
    return true;
  }
  return false;
}

void CglLKCI::generateCuts(const OsiSolverInterface& si, OsiCuts& cs, const CglTreeInfo info) {
  if (si.getNumCols() == 0 || si.getNumRows() == 0) return;
  const double* colSol = si.getColSolution();
  int numCols = si.getNumCols();
  int numRows = si.getNumRows();

  const char* colType = si.getColType(true);
  const double* colLB = si.getColLower();
  const double* colUB = si.getColUpper();

  const CoinPackedMatrix* matrixByRow = si.getMatrixByRow();
  const char* rowSense = si.getRowSense();
  const double* rowRHS = si.getRightHandSide();
  const double* rowRange = si.getRowRange();
  const double* rowActivity = si.getRowActivity();

  int numCoverRowsDetected = 0;
  CoinKnapsackRow knapsackRow(numCols, colType, colLB, colUB, 1e-7, si.getInfinity());

  for (int r = 0; r < numRows; ++r) {
    if (max_cover_rows_ >= 0 && numCoverRowsDetected >= max_cover_rows_) break;

    const CoinShallowPackedVector& row = matrixByRow->getVector(r);
    int rowNnz = row.getNumElements();
    if (max_row_nnz_ > 0 && rowNnz > max_row_nnz_) continue;

    const int* indices = row.getIndices();
    const double* elements = row.getElements();

    double multiplier[2];
    double adjustedRHS[2];
    int numIterations = CoinKnapsackRow::rowIterations(rowSense[r], rowRHS[r], rowRange[r], multiplier, adjustedRHS);

    for (int k = 0; k < numIterations; ++k) {
      if (rowActivity) {
        double slack = adjustedRHS[k] - multiplier[k] * rowActivity[r];
        if (slack >= 1.0 - 1e-7) {
          continue;
        }
      }

      knapsackRow.processRow(indices, elements, rowNnz, rowSense[r], multiplier[k], adjustedRHS[k]);
      if (knapsackRow.isUnbounded()) continue;

      size_t numTerms = knapsackRow.nzs();
      if (numTerms == 0) continue;

      const CoinTerm* terms = knapsackRow.columns();
      double b = knapsackRow.rhs();

      double sum_a = 0.0;
      for (size_t j = 0; j < numTerms; ++j) {
        sum_a += terms[j].value;
      }
      double d_double = sum_a - b;

      if (d_double <= 1e-9) continue;

      // Check if coefficients and rhs are near-integer
      bool ok = true;
      for (size_t j = 0; j < numTerms; ++j) {
        if (terms[j].value <= 1e-9 || std::abs(terms[j].value - std::round(terms[j].value)) > 1e-9) {
          ok = false;
          break;
        }
      }
      if (!ok || std::abs(d_double - std::round(d_double)) > 1e-9) continue;

      // Compute LP activity and number of fractional variables to prune early
      double lp_lhs = 0.0;
      int numFractional = 0;
      for (size_t j = 0; j < numTerms; ++j) {
        double val = (terms[j].index < numCols) ? (1.0 - colSol[terms[j].index]) : colSol[terms[j].index - numCols];
        double z_val = 1.0 - val;
        lp_lhs += terms[j].value * z_val;
        if (z_val > 1e-5 && z_val < 1.0 - 1e-5) {
          numFractional++;
        }
      }

      if (numFractional == 0 || lp_lhs <= b - 1.0 + 1e-7) {
        continue;
      }

      numCoverRowsDetected++;

      // Populate input vectors for separator
      std::vector<long long> a(numTerms);
      std::vector<double> row_x(numTerms);
      for (size_t j = 0; j < numTerms; ++j) {
        a[j] = (long long)std::round(terms[j].value);
        row_x[j] = (terms[j].index < numCols) ? (1.0 - colSol[terms[j].index]) : colSol[terms[j].index - numCols];
      }
      long long d = (long long)std::round(d_double);

      // Run separator
      LKCICandidate cand;
      if (separate_lkci(a, d, row_x, cand, r * 2 + k)) {
        // Build cut in terms of original variables x
        std::map<int, double> coeff;
        double cut_rhs_offset = 0.0;

        auto add_coeff_z = [&](int p, double coeff_z) {
          int index = terms[p].index;
          if (index < numCols) {
            // z_p = 1 - x_v, term is coeff_z * (1 - x_v) = coeff_z - coeff_z * x_v
            int col = index;
            coeff[col] -= coeff_z;
            cut_rhs_offset += coeff_z;
          } else {
            // z_p = x_v, term is coeff_z * x_v
            int col = index - numCols;
            coeff[col] += coeff_z;
          }
        };

        for (int p : cand.R) {
          add_coeff_z(p, (double)a[p]);
        }
        for (int p : cand.L) {
          add_coeff_z(p, (double)cand.D);
        }
        for (int p : cand.S) {
          add_coeff_z(p, cand.gamma[p]);
        }

        std::vector<int> cut_indices;
        std::vector<double> cut_values;
        for (const auto& kv : coeff) {
          if (std::abs(kv.second) > 1e-9) {
            cut_indices.push_back(kv.first);
            cut_values.push_back(kv.second);
          }
        }

        double sum_gamma_S = 0.0;
        for (int p : cand.S) {
          sum_gamma_S += cand.gamma[p];
        }
        double lower = (double)cand.D + sum_gamma_S - cut_rhs_offset;

        OsiRowCut rc;
        rc.setRow(cut_indices.size(), cut_indices.data(), cut_values.data());
        rc.setLb(lower);
        rc.setUb(si.getInfinity());
        cs.insertIfNotDuplicate(rc);
      }
    }
  }
}
