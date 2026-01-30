#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_permutations_cpp(
    IntegerVector vid_pairs_idx,    // indices into var_ids vector (0-based in C++)
    NumericVector a_pairs,          // annotation weights
    IntegerVector n_pairs,          // number of pairs per variable
    IntegerVector slot_count,       // slots per variable (n_pairs or 1)
    NumericVector base_weights,     // ranking weights by position
    int V,                          // total number of variables
    int n_perm,
    int seed
) {
  int n_pairs_total = vid_pairs_idx.size();
  
  if (n_pairs_total == 0 || V == 0) {
    return NumericVector(n_perm, NA_REAL);
  }
  
  // Precompute total_slots
  int total_slots = 0;
  for (int i = 0; i < V; i++) {
    total_slots += slot_count[i];
  }
  
  // Result vector
  NumericVector perm_ES(n_perm);
  
  // Set seed for reproducibility using R's set.seed
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
  
  // Permutation loop
  for (int b = 0; b < n_perm; b++) {
    // Shuffle variable order (1-indexed from R, convert to 0-indexed)
    IntegerVector perm_order = Rcpp::sample(V, V, false);
    for (int i = 0; i < V; i++) perm_order[i] -= 1;
    
    // Build permuted skeleton
    std::vector<int> n_pairs_perm(V);
    std::vector<int> slot_count_perm(V);
    for (int i = 0; i < V; i++) {
      n_pairs_perm[i] = n_pairs[perm_order[i]];
      slot_count_perm[i] = slot_count[perm_order[i]];
    }
    
    // Rebuild pair_positions and miss_positions for this permutation
    std::vector<int> pair_pos_perm;
    std::vector<int> miss_pos_perm;
    int current_pos = 0;
    for (int i = 0; i < V; i++) {
      if (n_pairs_perm[i] > 0) {
        for (int j = 0; j < n_pairs_perm[i]; j++) {
          pair_pos_perm.push_back(current_pos + j);
        }
        current_pos += n_pairs_perm[i];
      } else {
        miss_pos_perm.push_back(current_pos);
        current_pos += 1;
      }
    }
    
    // Create permuted var2pos mapping
    std::vector<int> var2pos_perm(V);
    for (int i = 0; i < V; i++) {
      var2pos_perm[perm_order[i]] = i;
    }
    
    // Compute pair weights using permuted mapping
    std::vector<double> pw;
    std::vector<int> valid_idx;
    for (int i = 0; i < n_pairs_total; i++) {
      int var_idx = vid_pairs_idx[i];
      int pos = var2pos_perm[var_idx];
      double weight = base_weights[pos];
      double pair_weight = a_pairs[i] * weight;
      
      if (!NumericVector::is_na(pair_weight) && pair_weight > 0) {
        pw.push_back(pair_weight);
        valid_idx.push_back(i);
      }
    }
    
    int n_valid = pw.size();
    if (n_valid == 0) {
      perm_ES[b] = NA_REAL;
      continue;
    }
    
    // Sort pairs by weight (descending)
    std::vector<int> ord(n_valid);
    for (int i = 0; i < n_valid; i++) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
      return pw[a] > pw[b];
    });
    
    int n_fill = std::min((int)pair_pos_perm.size(), n_valid);
    if (n_fill == 0) {
      perm_ES[b] = 0.0;
      continue;
    }
    
    // Compute hit total
    double hit_total = 0.0;
    for (int i = 0; i < n_fill; i++) {
      hit_total += pw[ord[i]];
    }
    
    if (hit_total <= 0 || !std::isfinite(hit_total)) {
      perm_ES[b] = NA_REAL;
      continue;
    }
    
    // Build step vector
    std::vector<double> step(total_slots, 0.0);
    for (int i = 0; i < n_fill; i++) {
      step[pair_pos_perm[i]] = pw[ord[i]] / hit_total;
    }
    
    if (miss_pos_perm.size() > 0) {
      double miss_step = -1.0 / miss_pos_perm.size();
      for (size_t i = 0; i < miss_pos_perm.size(); i++) {
        step[miss_pos_perm[i]] = miss_step;
      }
    }
    
    // Cumulative sum and max ES
    double cumsum = 0.0;
    double max_es = 0.0;
    for (int i = 0; i < total_slots; i++) {
      cumsum += step[i];
      if (cumsum > max_es) max_es = cumsum;
    }
    
    perm_ES[b] = max_es;
  }
  
  return perm_ES;
}