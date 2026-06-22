#ifndef CglLKCI_H
#define CglLKCI_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include "CglCutGenerator.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

/**
 * Lifting Knapsack Cover Inequalities (LKCI) Heuristic Cut Generator
 */
class CGLLIB_EXPORT CglLKCI : public CglCutGenerator {
public:
  /**@name Generate Cuts */
  //@{
  /** Generate LKCI cuts for the model accessed through the solver interface.
      Insert generated cuts into the cut set cs.
  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
                            const CglTreeInfo info = CglTreeInfo());
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor
  CglLKCI();

  /// Copy constructor
  CglLKCI(const CglLKCI &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator
  CglLKCI & operator=(const CglLKCI & rhs);

  /// Destructor
  virtual ~CglLKCI();

  /// Create C++ lines to get to current state
  virtual std::string generateCpp(FILE * fp);
  //@}

  /**@name Gets and Sets */
  //@{
  inline int getBeamWidth() const { return beam_width_; }
  inline void setBeamWidth(int value) { beam_width_ = value; }

  inline int getRandomStarts() const { return random_starts_; }
  inline void setRandomStarts(int value) { random_starts_ = value; }

  inline int getLocalRounds() const { return local_rounds_; }
  inline void setLocalRounds(int value) { local_rounds_ = value; }

  inline int getMovePool() const { return move_pool_; }
  inline void setMovePool(int value) { move_pool_ = value; }

  inline int getSeed() const { return seed_; }
  inline void setSeed(int value) { seed_ = value; }

  inline double getTolerance() const { return tol_; }
  inline void setTolerance(double value) { tol_ = value; }

  inline int getMaxCoverRows() const { return max_cover_rows_; }
  inline void setMaxCoverRows(int value) { max_cover_rows_ = value; }

  inline int getMaxRowNnz() const { return max_row_nnz_; }
  inline void setMaxRowNnz(int value) { max_row_nnz_ = value; }
  //@}

private:
  // Heuristic parameters
  int beam_width_;
  int random_starts_;
  int local_rounds_;
  int move_pool_;
  int seed_;
  double tol_;
  int max_cover_rows_;
  int max_row_nnz_;

  // Helper structures for candidate evaluation
  struct LKCICandidate {
    std::vector<int> S;
    long long D;
    double delta;
    std::vector<int> L;
    std::vector<int> R;
    std::vector<double> gamma;
    double lhs;
    double rhs;
  };

  struct BeamEntry {
    double delta;
    std::vector<bool> mask;
    LKCICandidate cand;

    bool operator<(const BeamEntry& other) const {
      return delta > other.delta; // Descending order
    }
  };

  // Helper methods
  static double h_value(long long r, long long D, const std::vector<long long>& T, const std::vector<long long>& prefix_min);
  static std::vector<double> h_values(const std::vector<long long>& r, long long D, const std::vector<long long>& b_desc);
  
  static bool evaluate_lkci(
      const std::vector<long long>& a,
      long long d,
      const std::vector<double>& x,
      const std::vector<bool>& S_mask,
      LKCICandidate& cand,
      double tol,
      bool require_positive);

  static std::vector<int> topk_indices(const std::vector<double>& scores, const std::vector<bool>& mask, int k);

  static void add_prefix_candidates(
      const std::vector<long long>& a,
      long long d,
      const std::vector<double>& x,
      const std::vector<int>& order,
      const std::vector<long long>& targets,
      std::vector<BeamEntry>& beam,
      int beam_width,
      double tol);

  bool separate_lkci(
      const std::vector<long long>& a,
      long long d,
      const std::vector<double>& x,
      LKCICandidate& cand_out,
      int seed_offset);

  static bool run_exact_best_move(
      const std::vector<long long>& a,
      long long d,
      const std::vector<double>& x,
      const std::vector<long long>& targets_vec,
      std::vector<bool>& mask,
      LKCICandidate& cand,
      double tol,
      int move_pool,
      double& best_gain);
};

#endif
