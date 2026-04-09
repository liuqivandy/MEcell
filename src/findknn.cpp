#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>
#include <numeric>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// ---------------------------------------------------------------------------
// KD-tree with OpenMP-parallel queries (replacement for nabor::knn)
// ---------------------------------------------------------------------------
// Data layout: R column-major matrix (n x d).  Point i, dimension j is at
// data[j * n + i].  The tree is built single-threaded; each query is
// independent, so queries are parallelised with OpenMP.

static const int KD_LEAF_SIZE = 32;

struct KDNode
{
  int split_dim; // -1 for leaf
  double split_val;
  int left, right; // child indices (-1 = none)
  int leaf_start;  // start offset into reordered index array
  int leaf_count;  // number of points in this leaf
};

class KDTree
{
public:
  std::vector<KDNode> nodes;
  std::vector<int> idx; // point indices, reordered during build
  const double *data;   // pointer to column-major R matrix
  int n, d;

  KDTree(const double *data, int n, int d) : data(data), n(n), d(d)
  {
    idx.resize(n);
    std::iota(idx.begin(), idx.end(), 0);
    nodes.reserve(4 * n / KD_LEAF_SIZE + 16);
    build(0, n);
  }

  // Build sub-tree over idx[start .. start+count-1].  Returns node index.
  int build(int start, int count)
  {
    int nid = (int)nodes.size();
    nodes.push_back(KDNode());

    if (count <= KD_LEAF_SIZE)
    {
      nodes[nid].split_dim = -1;
      nodes[nid].leaf_start = start;
      nodes[nid].leaf_count = count;
      nodes[nid].left = nodes[nid].right = -1;
      return nid;
    }

    // Pick dimension with largest spread
    int best_dim = 0;
    double max_spread = 0.0;
    for (int dim = 0; dim < d; dim++)
    {
      double mn = std::numeric_limits<double>::max();
      double mx = std::numeric_limits<double>::lowest();
      const double *col = data + (size_t)dim * n;
      for (int i = 0; i < count; i++)
      {
        double v = col[idx[start + i]];
        if (v < mn)
          mn = v;
        if (v > mx)
          mx = v;
      }
      double spread = mx - mn;
      if (spread > max_spread)
      {
        max_spread = spread;
        best_dim = dim;
      }
    }

    // Partition around median using nth_element
    int mid = count / 2;
    const double *split_col = data + (size_t)best_dim * n;
    std::nth_element(idx.begin() + start,
                     idx.begin() + start + mid,
                     idx.begin() + start + count,
                     [split_col](int a, int b)
                     {
                       return split_col[a] < split_col[b];
                     });

    nodes[nid].split_dim = best_dim;
    nodes[nid].split_val = split_col[idx[start + mid]];
    nodes[nid].leaf_start = -1;
    nodes[nid].leaf_count = 0;

    int left = build(start, mid);
    int right = build(start + mid, count - mid);
    nodes[nid].left = left;
    nodes[nid].right = right;
    return nid;
  }

  // KNN query: maintains a max-heap of size k (largest dist on top).
  // Called recursively; prunes branches whose splitting plane is farther
  // than the current k-th distance.
  void query(const double *q, int k,
             std::priority_queue<std::pair<double, int>> &heap,
             int nid) const
  {
    const KDNode &nd = nodes[nid];

    if (nd.split_dim == -1)
    {
      // Leaf – brute force against all points in this bucket
      for (int i = 0; i < nd.leaf_count; i++)
      {
        int pi = idx[nd.leaf_start + i];
        double dist2 = 0.0;
        for (int dim = 0; dim < d; dim++)
        {
          double diff = q[dim] - data[(size_t)dim * n + pi];
          dist2 += diff * diff;
        }
        if ((int)heap.size() < k)
        {
          heap.push({dist2, pi});
        }
        else if (dist2 < heap.top().first)
        {
          heap.pop();
          heap.push({dist2, pi});
        }
      }
      return;
    }

    double diff = q[nd.split_dim] - nd.split_val;
    int first = (diff <= 0.0) ? nd.left : nd.right;
    int second = (diff <= 0.0) ? nd.right : nd.left;

    query(q, k, heap, first);

    double plane_dist2 = diff * diff;
    if ((int)heap.size() < k || plane_dist2 < heap.top().first)
    {
      query(q, k, heap, second);
    }
  }
};

// [[Rcpp::export]]
List knn_rcpp(NumericMatrix data, int k, int nthreads = 0)
{
  int n = data.nrow();
  int d = data.ncol();
  const double *data_ptr = REAL(data);

  // Build kd-tree (single-threaded, O(n log n))
  KDTree tree(data_ptr, n, d);

  // Allocate output (column-major: n x k)
  IntegerMatrix nn_idx(n, k);
  NumericMatrix nn_dists(n, k);
  int *idx_out = INTEGER(nn_idx);
  double *dist_out = REAL(nn_dists);

// Parallel queries – each point's search is independent
#ifdef _OPENMP
  if (nthreads > 0)
    omp_set_num_threads(nthreads);
#pragma omp parallel for schedule(dynamic, 64)
#endif
  for (int i = 0; i < n; i++)
  {
    // Gather query coordinates (strided in column-major)
    std::vector<double> q(d);
    for (int dim = 0; dim < d; dim++)
    {
      q[dim] = data_ptr[(size_t)dim * n + i];
    }

    std::priority_queue<std::pair<double, int>> heap;
    tree.query(q.data(), k, heap, 0);

    // Heap gives largest-first; pop into output in reverse order
    int found = (int)heap.size();
    for (int j = found - 1; j >= 0; j--)
    {
      dist_out[i + (size_t)j * n] = std::sqrt(heap.top().first);
      idx_out[i + (size_t)j * n] = heap.top().second + 1; // 1-based
      heap.pop();
    }
  }

  return List::create(Named("nn.idx") = nn_idx,
                      Named("nn.dists") = nn_dists);
}

// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List find_knn_rcpp(NumericMatrix mat, IntegerMatrix candidates_mat)
{
  int n_samples = mat.ncol();
  int n_genes = mat.nrow();
  int n_cand = candidates_mat.ncol();

  // Output matrices — allocated before parallel region
  NumericMatrix all_distances(n_samples, n_cand);
  IntegerMatrix all_indices(n_samples, n_cand);

  // Raw pointers for thread-safe, zero-overhead access
  const double *mat_ptr = REAL(mat);
  const int *cand_ptr = INTEGER(candidates_mat);
  double *dist_ptr = REAL(all_distances);
  int *idx_ptr = INTEGER(all_indices);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
  for (int i = 0; i < n_samples; i++)
  {
    // Direct pointer to query column (column-major: column i starts at i*n_genes)
    const double *query = mat_ptr + (size_t)i * n_genes;

    std::vector<std::pair<double, int>> dist_idx(n_cand);

    for (int j = 0; j < n_cand; j++)
    {
      // candidates_mat is n_samples x n_cand, column-major
      int c_idx = cand_ptr[i + (size_t)j * n_samples] - 1; // 0-based
      const double *cand_col = mat_ptr + (size_t)c_idx * n_genes;

      // Compute squared Euclidean distance (skip sqrt for sorting)
      double dist = 0.0;
      for (int g = 0; g < n_genes; g++)
      {
        double diff = cand_col[g] - query[g];
        dist += diff * diff;
      }

      dist_idx[j] = std::make_pair(dist, c_idx + 1); // 1-based for R
    }

    // Sort by squared distance (monotonic with Euclidean distance)
    std::sort(dist_idx.begin(), dist_idx.end(),
              [](const std::pair<double, int> &a, const std::pair<double, int> &b)
              {
                return a.first < b.first;
              });

    // Write results — compute sqrt only here for final output
    for (int j = 0; j < n_cand; j++)
    {
      dist_ptr[i + (size_t)j * n_samples] = std::sqrt(dist_idx[j].first);
      idx_ptr[i + (size_t)j * n_samples] = dist_idx[j].second;
    }
  }

  return List::create(Named("distances") = all_distances,
                      Named("indices") = all_indices);
}

// C++ implementation of CalMEI to replace the R for-loop
// Takes the sparse dgCMatrix components directly
// [[Rcpp::export]]
NumericVector cal_mei_rcpp(IntegerVector p, IntegerVector i_vec, NumericVector x,
                           IntegerVector cellcluster, int topn)
{
  int ncell = p.size() - 1;
  NumericVector mei(ncell);

  for (int cell = 0; cell < ncell; cell++)
  {
    int start = p[cell];
    int end = p[cell + 1];
    int n_neighbors = end - start;
    if (n_neighbors == 0)
    {
      mei[cell] = 1.0; // no neighbors = max MEI
      continue;
    }

    // Collect neighbor weights and indices
    int top_k = std::min(topn, n_neighbors);
    std::vector<std::pair<double, int>> weight_idx(n_neighbors);
    for (int k = 0; k < n_neighbors; k++)
    {
      weight_idx[k] = std::make_pair(x[start + k], i_vec[start + k]);
    }

    // Partial sort: only need the top-k by descending weight
    std::partial_sort(weight_idx.begin(), weight_idx.begin() + top_k, weight_idx.end(),
                      [](const std::pair<double, int> &a, const std::pair<double, int> &b)
                      {
                        return a.first > b.first;
                      });

    // Count how many of the top-k share the cluster of the top neighbor
    int ref_cluster = cellcluster[weight_idx[0].second]; // 0-based index into cellcluster
    int support = 0;
    for (int k = 0; k < top_k; k++)
    {
      if (cellcluster[weight_idx[k].second] == ref_cluster)
      {
        support++;
      }
    }

    mei[cell] = (double)(topn - support) / topn;
  }

  return mei;
}
