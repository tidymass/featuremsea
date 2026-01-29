#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_permutations_cpp(
    IntegerVector vid_pairs_idx,     // Feature indices (0-based in C++)
    IntegerVector met_ids_encoded,   // Metabolite IDs (encoded as integers)
    NumericVector a_pairs,           // Annotation scores (after feature-level dedup)
    NumericVector base_weights,      // Ranking weights by position
    int V,                           // Total number of variables (features)
    int n_perm,                      // Number of permutations
    int seed                         // Random seed
) {
  int n_pairs_total = vid_pairs_idx.size();
  
  if (n_pairs_total == 0 || V == 0) {
    return NumericVector(n_perm, NA_REAL);
  }
  
  // Result vector
  NumericVector perm_ES(n_perm);
  
  // Set seed for reproducibility
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
  
  // Permutation loop
  for (int b = 0; b < n_perm; b++) {
    // 1. Shuffle variable order (sample returns 1-indexed, convert to 0-indexed)
    IntegerVector perm_order = Rcpp::sample(V, V, false);
    for (int i = 0; i < V; i++) perm_order[i] -= 1;
    
    // 2. Build permuted var2pos mapping
    std::vector<int> var2pos_perm(V);
    for (int i = 0; i < V; i++) {
      var2pos_perm[perm_order[i]] = i;
    }
    
    // 3. Calculate new rank positions for all pairs
    std::vector<int> rank_positions(n_pairs_total);
    for (int i = 0; i < n_pairs_total; i++) {
      rank_positions[i] = var2pos_perm[vid_pairs_idx[i]];
    }
    
    // ============================================================
    // 4. METABOLITE-LEVEL DEDUPLICATION
    // Keep only the best-ranked feature for each metabolite
    // ============================================================
    std::unordered_map<int, int> best_feature_per_met;  // met_id -> index in pairs
    
    for (int i = 0; i < n_pairs_total; i++) {
      int met_id = met_ids_encoded[i];
      int curr_rank = rank_positions[i];
      
      // If this metabolite hasn't been seen, or current feature has better rank
      if (best_feature_per_met.find(met_id) == best_feature_per_met.end()) {
        best_feature_per_met[met_id] = i;
      } else {
        int best_idx = best_feature_per_met[met_id];
        if (curr_rank < rank_positions[best_idx]) {  // Lower rank = better
          best_feature_per_met[met_id] = i;
        }
      }
    }
    
    // 5. Calculate pair weights for deduplicated pairs only
    std::vector<double> pw;
    std::vector<int> valid_pos;  // Positions in the ranking
    
    for (auto& kv : best_feature_per_met) {
      int idx = kv.second;
      int pos = rank_positions[idx];
      double weight = base_weights[pos];
      double pair_weight = a_pairs[idx] * weight;
      
      if (!NumericVector::is_na(pair_weight) && pair_weight > 0) {
        pw.push_back(pair_weight);
        valid_pos.push_back(pos);
      }
    }
    
    int n_valid = pw.size();
    if (n_valid == 0) {
      perm_ES[b] = NA_REAL;
      continue;
    }
    
    // 6. Sort pairs by weight (descending)
    std::vector<int> ord(n_valid);
    for (int i = 0; i < n_valid; i++) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
      return pw[a] > pw[b];
    });
    
    // 7. Calculate hit total
    double hit_total = 0.0;
    for (int i = 0; i < n_valid; i++) {
      hit_total += pw[ord[i]];
    }
    
    if (hit_total <= 0 || !std::isfinite(hit_total)) {
      perm_ES[b] = NA_REAL;
      continue;
    }
    
    // 8. Build step vector for ES calculation
    // Total slots = V (one slot per feature in the ranking)
    std::vector<double> step(V, 0.0);
    
    // Mark which positions are hits
    std::vector<bool> is_hit(V, false);
    for (int i = 0; i < n_valid; i++) {
      int pos = valid_pos[ord[i]];
      step[pos] = pw[ord[i]] / hit_total;
      is_hit[pos] = true;
    }
    
    // Calculate number of misses
    int n_miss = 0;
    for (int i = 0; i < V; i++) {
      if (!is_hit[i]) n_miss++;
    }
    
    // Assign miss penalty
    if (n_miss > 0) {
      double miss_step = -1.0 / n_miss;
      for (int i = 0; i < V; i++) {
        if (!is_hit[i]) {
          step[i] = miss_step;
        }
      }
    }
    
    // 9. Cumulative sum and find max ES
    double cumsum = 0.0;
    double max_es = 0.0;
    for (int i = 0; i < V; i++) {
      cumsum += step[i];
      if (cumsum > max_es) max_es = cumsum;
    }
    
    perm_ES[b] = max_es;
  }
  
  return perm_ES;
}