#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
List find_knn_rcpp(NumericMatrix mat, IntegerMatrix candidates_mat) {
  int n_samples = mat.ncol();
  int n_genes = mat.nrow();
  int n_cand = candidates_mat.ncol();

  // Output matrices
  NumericMatrix all_distances(n_samples, n_cand);
  IntegerMatrix all_indices(n_samples, n_cand);

  for (int i = 0; i < n_samples; i++) {
    NumericVector query = mat(_, i);

    std::vector<std::pair<double, int>> dist_idx(n_cand);

    for (int j = 0; j < n_cand; j++) {
      int c_idx = candidates_mat(i, j) - 1;  // Convert to 0-based
      double dist = 0.0;

      for (int g = 0; g < n_genes; g++) {
        double diff = mat(g, c_idx) - query[g];
        dist += diff * diff;
      }

      dist_idx[j] = std::make_pair(std::sqrt(dist), c_idx + 1); // back to 1-based
    }

    // Sort candidates by distance
    std::sort(dist_idx.begin(), dist_idx.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                return a.first < b.first;
              });

    // Fill output matrices
    for (int j = 0; j < n_cand; j++) {
      all_distances(i, j) = dist_idx[j].first;
      all_indices(i, j) = dist_idx[j].second;
    }
  }

  return List::create(Named("distances") = all_distances,
                      Named("indices") = all_indices);
}
