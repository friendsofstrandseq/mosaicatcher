#include <cassert>
#include <iostream>
#include <vector>
#include <limits>
#include <math.h>  // M_PI
#include <iomanip> // setw
#include <numeric> // partial_sum

/**
 * @file
 * @defgroup segmentation Segmentation algorithm
 * Summary of how segments across binned Strand-seq data of multiple cells are
 * found.
 *
 * ### Overview of the segmentation algorithm
 *
 * We apply an optimal multivariate segmentation algorithm to find segments
 * across all cells simultaneously.
 * @todo Write description
 */


template<class T>
using Matrix = std::vector<std::vector<T> >;

/**
 * @fn template <typename TMat> void print_mat(TMat const & G)
 * @ingroup segmentation
 * Print a Matrix (vector<vector>).
 *
 * @param G 2D matrix[row_index][column_index]
 */
template <typename TMat> void print_mat(TMat const & G) {
    for (auto row : G) {
        for (auto x : row)
            std::cout << std::setw(10) << x;
        std::cout << std::endl;
    }
}

/**
 * @fn bool segmentation(Matrix<float> const & cost, int max_cp, Matrix<int> & breakpoints, std::vector<float> & log_lik)
 * @ingroup segmentation
 * Find optimal segmentation based on a cost matrix.
 *  
 * This is a dynamic programming algorithm to find an optimal segmentation
 * in a cost matrix (see calculate_cost_matrix). It calculates a max_cp x N
 * matrix of optimal segmentation cost (internally) and breakpoints (returned
 * via `breakpoints`).
 *
 * ### Optimal cost matrix `dp` of size N x max_cp
 *
 * `dp[n, cp]` stores the cost of the optimal segmentation with
 * `cp` breakpoints of segment `[0, n]`. This is calculated via
 * dynamic programming from the optimal cost of a segmentation
 * of `[0, k-1]` with `cp-1` breakpoints (note that 0 <= k < n) plus
 * the actual cost of segment `[k, n]` from the "cost" matrix.
 *
 * @param cost Cost matrix (see `calculate_cost_matrix`)
 * @param max_cp Maximum number of breakpoints (<= N)
 * @param breakpoints Breakpoints for k=1..max_cp will be written in here
 * @param log_lik Something similar to BIC is written in here. Might be useful for model selection.
 */
bool segmentation(Matrix<float> const & cost,
                  int max_cp,
                  Matrix<int> & breakpoints,
                  std::vector<float> & log_lik)
{
    // Determine max_k and N from cost matrix
    unsigned max_k  = (unsigned)cost.size();
    if (max_k < 5 || max_k > 5000) {
        std::cerr << "[ERROR] unreasonable value for max_k: " << max_k << std::endl;
        return false;
    }
    unsigned N      = (unsigned)cost[0].size();
    if (N < 5 || N > 25000) {
        std::cerr << "[ERROR] unreasonable value for N: " << N << std::endl;
        return false;
    }
    if (max_cp > N || max_k > N) {
        std::cerr << "[ERROR] max_cp or max_k cannot be larger than N: N=" << N << ", max_cp=" << max_cp << ", max_k=" << max_k << std::endl;
        return false;
    }

    // cost of optimal segmentation matrix
    Matrix<float> dp(N, std::vector<float>(max_cp, -1));

    // backtrack matrix mt
    Matrix<int>    mt(N, std::vector<int>(max_cp - 1, -1));

    // Initialization: dp[n][0] is just cost[n][0]
    for (unsigned r = 0; r < max_k; ++r)
        dp[r][0] = cost[r][0];
    for (unsigned r = max_k; r < N; ++r)
        dp[r][0] = -1;


    // Dynamic programming: column-wise
    for (unsigned cp=1; cp < max_cp; ++cp)
    {
        for (unsigned j=0; j<N; ++j)
        {
            float   z_min = -1;
            unsigned i_min =  j;
            int      j0    = (j<max_k) ? j : max_k;

            // Iterate over possible new mt between 0 and j
            for (unsigned k=0; k<j0; ++k)
            {
                // Best segmentation from 0 to j-k-1 with cp-1 change points
                float mI_prev_col = dp[j-k-1][cp-1];

                // Cost of segment from j-k to j
                float cost_seg    = cost[k][j-k];

                // Keep value with minimal cost and set mt matriv=x
                // Note that values < 0 represent positive infinity
                if (mI_prev_col > 0 && cost_seg > 0)
                {
                    if (z_min < 0 || mI_prev_col + cost_seg < z_min) {
                        z_min = mI_prev_col + cost_seg;
                        i_min = j-k;
                    }
                }
            } /* for k */
            dp[j][cp]   = z_min;
            mt[j][cp-1] = i_min;
        } /* for j */
    } /* for cp */

    // breakpoints: Row cp contains the breakpoints for a segmentation with cp breakpoints.
    // breakpoints[cp][cp] is always N
    log_lik     = std::vector<float>(max_cp);
    breakpoints = Matrix<int>(max_cp, std::vector<int>(max_cp, -1)); // important: initialize to -1

    for (unsigned cp = 0; cp < max_cp; ++cp)
    {
        /* Calculate J, the log-likelihood. */
        float z = dp[N-1][cp];
        log_lik[cp] = (z > 0) ? -(float)N / 2.0 * (1 + log(2*M_PI) + log(z / N)) : -1e10;  // 1e10 as a really really low value

        // Backtrack to get breakpoints
        // i is always the oosition of the changepoint to the right (???)
        int i = (int)N;
        breakpoints[cp][cp] = N;
        for (int j = cp - 1; j >= 0; --j) {
            i = mt[i-1][j];
            breakpoints[cp][j] = i;
        }
    }
    std::cout << "Vector logLik" << std::endl;
    for (auto x : log_lik)
        std::cout << x << "\t";
    std::cout << std::endl;

    std::cout << "Matrix breakpoints" << std::endl;
    print_mat(breakpoints);

    return true;
}


/**
 * @fn Matrix<float> calculate_cost_matrix(std::vector<float> const & data, int max_k)
 * @ingroup segmentation
 * Calculate segmentation cost matrix on a single vector.
 *
 * Function re-implemented from Wolfgang Hubers tilingArray package
 * (https://github.com/Bioconductor-mirror/tilingArray/blob/master/R/costMatrix.R).
 * This is the version that runs on a single vector.
 *
 * @todo Profile this algorith, to see where time is spent.
 * @todo Potentially replace calculate_cost_matrix_single by a faster version using
 *       boost UBLAS vectors and matrices.
 * 
 * @param data Single data vector.
 * @param max_k Maximum number of bins for each segment. max_k <= N.
 */
Matrix<float> calculate_cost_matrix(std::vector<float> const & data,
                                                int max_k)
{
    // See https://www.bioconductor.org/packages/devel/bioc/vignettes/tilingArray/inst/doc/costMatrix.pdf
    // for an explanation of the algebra
    unsigned N = (unsigned)data.size();
    std::vector<float> cr(N);
    std::vector<float> cq(N);

    // cr = cumsum(data)
    std::partial_sum(data.begin(), data.end(), cr.begin());                             // O(N)

    // cq = cumsum(data^2)
    std::transform(data.begin(), data.end(), cq.begin(), [](float u){return u*u;});     // O(N)
    std::partial_sum(cq.begin(), cq.end(), cq.begin());                                 // O(N)

    // Cost matrix G
    Matrix<float> G(max_k, std::vector<float>(N, -1));                                  // O(N * maxk)

    // Initialization:
    for (unsigned k = 0; k < max_k; ++k)                                                // O(maxk)
        G[k][0] = cq[k] - cr[k]*cr[k]/(float)(k+1);

    // Iterate per row (k)
    for (unsigned k = 1; k <= max_k && k <= N-1; ++k)                                   // O( maxk * ...
    {
        std::vector<float> cqk(N-k, -1);                                                //           ... N
        std::vector<float> crk(N-k, -1);                                                //           ... N

        for (unsigned i = 0; i < N-k; ++i) {                                            //           ... N
            cqk[i] = cq[i+k] - cq[i];
            crk[i] = cr[i+k] - cr[i];
        }

        // Iterate columns j = 1...N-k
        for (unsigned j = 1; j < N-k+1; ++j) {                                          //           ... N
            unsigned i = j - 1;
            assert(cqk[i]>=0);
            assert(crk[i]>=0);
            G[k-1][j] = cqk[i] - crk[i]*crk[i]/(float)k;
        }
    }

    return G;                                                                           // TOTAL: O(N * maxk) <= O(N^2)
}


/**
 * @fn calculate_cost_matrix(Matrix<float> const & data, int max_k, Matrix<float> & G)
 * @ingroup segmentation
 * Calculate segmentation cost matrix.
 *
 * Function re-implemented from Wolfgang Hubers tilingArray
 * (https://github.com/Bioconductor-mirror/tilingArray/blob/master/R/costMatrix.R).
 * The difference to the tilingArray package is that, unlike there,
 * a separate piecewise-constant function is fitted to every sample
 * (in the original version all samples are averaged into one).
 * This is done by calculating the cost matrix seperately on each sample
 * and then averaging across them. Complexity is quite high now, with
 * O(N * maxk * J) where J is the number of strands or samples.
 *
 * @todo Test for correctness (e.g. write unit test)
 *
 * @param data Matrix of data values, each row representing a single sample 
 *             (or strand). Note that all must have the same length and should
 *             be normalized beforehands!
 * @param max_k Number of bins of longest possible segment. This reduces the
 *              complexity from O(N^2) to O(N*max_k).
 * @param G Cost matrix. Reference will be resized and filled.
 * 
 * @addtogroup segmentation
 */
bool calculate_cost_matrix(Matrix<float> const & data,                         // TOTAL O ( J * N * maxk)
                                    int max_k,
                                    Matrix<float> & G)
{

    unsigned J  = (unsigned)data.size();        // number of cells or strands
    if (J < 1 || J > 2000) {
        std::cerr << "[ERROR] more than 2000 strands might be too much: " << J << std::endl;
        return false;
    }
    unsigned N      = (unsigned)data[0].size(); // length of signal
    if (N < 5 || N > 25000) {
        std::cerr << "[ERROR] unreasonable value for N: " << N << std::endl;
        return false;
    }
    if (max_k > N) {
        std::cerr << "[ERROR] max_k cannot be larger than N: N=" << N << ", max_k=" << max_k << std::endl;
        return false;
    }


    G = Matrix<float>(max_k, std::vector<float>(N, -1));
    for (unsigned s = 0; s < J; ++s) {

        // Get cost of this one sample
        Matrix<float> cost_j = calculate_cost_matrix(data[s], max_k);

        // Add to total cost
        for (unsigned i=0; i<max_k; ++i)
            for (unsigned j=0; j<N; ++j)
                G[i][j] += cost_j[i][j];

    }

    // Finally, divide by J
    for (unsigned i=0; i<max_k; ++i)
        for (unsigned j=0; j<N; ++j)
            G[i][j] /= (float)J;

    return true;
}



int main() {

    //                            1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
    //                          102  102  100  100  100  100  100  100  100  104  104  100  100  100  100
    //                          =========-----------------------------------==========-------------------
    Matrix<float> data   = {{100, 100, 100, 100, 100,  50,  50,  50,  50,  50, 100, 100, 100, 100, 100},
                            {  2,   2,   0,   0,   0,  50,  50,  50,  50,  54,   4,   0,   0,   0,   0}};

    // real values:
    unsigned max_cp = 4;
    unsigned max_k  = 15;

    // New Cost matrix
    Matrix<float> new_cost;
    calculate_cost_matrix(data, max_k, new_cost);
    std::cout << "NEW cost matrix (maxk x N):" << std::endl;
    print_mat(new_cost);

    // Find optimal segmentation
    Matrix<int> breakpoints;
    std::vector<float> log_lik;
    if (!segmentation(new_cost, max_cp, breakpoints, log_lik))
        std::cerr << "[ERROR] Segmentation failed" << std::endl;



} // main



