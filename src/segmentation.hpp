/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#include <cassert>
#include <iostream>
#include <vector>
#include <limits>
#include <math.h>  // M_PI
#include <iomanip> // setw
#include <numeric> // partial_sum
#include <chrono>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>

#include "version.hpp"
#include "utils.hpp"
#include "iocounts.hpp"


using interval::Interval;
using count::TGenomeCounts;
using count::Counter;

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
 *
 * See `main_segment` for a description of which steps are done.
 *
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

template <typename Printable>
std::ostream& operator<< (std::ostream& stream, Matrix<Printable> const & m) {
    if (m.size() == 0 || m[0].size() == 0)
        return (stream << "[EMPTY MATRIX]" << std::endl);
    for (unsigned r = 0; r < m.size(); ++r) {
        //stream << std::setw(8) << std::setprecision(2) << m[r][0];
        for (unsigned c = 0; c < m[r].size(); ++c)
            stream << std::setw(10) << std::setprecision(4) << m[r][c];
        stream << std::endl;
    }
    return stream;
}



/**
 * @fn bool optimal_segment_dp(Matrix<double> const & cost, int max_cp, Matrix<int> & breakpoints, std::vector<double> & sse)
 * @ingroup segmentation
 * Find optimal segmentation based on a cost matrix.
 *  
 * This is a dynamic programming algorithm to find an optimal segmentation
 * in a cost matrix (see `calculate_cost_matrix`). It calculates a max_cp x N
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
 * @param sse sum of squared error value written here.
 */
bool optimal_segment_dp(Matrix<double> const & cost,
                  int max_cp,
                  Matrix<int> & breakpoints,
                  std::vector<double> & sse)
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
    Matrix<double> dp(N, std::vector<double>(max_cp, -1));

    // backtrack matrix mt
    Matrix<int>    mt(N, std::vector<int>(max_cp - 1, -1));

    // Initialization: dp[i][0] is just cost[i][0]
    for (unsigned r = 0; r < max_k; ++r)
        dp[r][0] = cost[r][0];
    for (unsigned r = max_k; r < N; ++r)
        dp[r][0] = -1;


    // Dynamic programming: column-wise (number of change points).
    for (unsigned cp=1; cp < max_cp; ++cp)
    {
        // row j = position j
        for (unsigned j=0; j<N; ++j)
        {
            double   z_min = -1;
            unsigned i_min =  j;
            int      j0    = (j<max_k) ? j : max_k;

            // Iterate over possible new mt between 0 and j
            for (unsigned k=0; k<j0; ++k)
            {
                // Best segmentation from 0 to j-k-1 with cp-1 change points
                double mI_prev_col = dp[j-k-1][cp-1];

                // Cost of segment from j-k to j
                double cost_seg    = cost[k][j-k];

                // Keep value with minimal cost and set mt matriv=x
                // Note that values < 0 represent positive infinity
                if (mI_prev_col >= 0 && cost_seg >= 0)
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
    sse     = std::vector<double>(max_cp);
    breakpoints = Matrix<int>(max_cp, std::vector<int>(max_cp, -1)); // important: initialize to -1

    for (unsigned cp = 0; cp < max_cp; ++cp)
    {
        // just write down SSE
        double z = dp[N-1][cp];
        //sse[cp] = (z > 0) ? -(double)N / 2.0 * (1 + log(2*M_PI) + log(z / N)) : -1e10;  // 1e10 as a really really low value
        sse[cp] = (z >= 0) ? z : 1e10;  // 1e10 as a really really high value

        // Backtrack to get breakpoints
        // i is always the oosition of the changepoint to the right (???)
        int i = (int)N;
        breakpoints[cp][cp] = N;
        for (int j = cp - 1; j >= 0; --j) {
            i = mt[i-1][j];
            breakpoints[cp][j] = i;
        }
    }
    return true;
}


/**
 * @fn Matrix<double> calculate_cost_matrix(std::vector<double> const & data, int max_k)
 * @ingroup segmentation
 * Calculate segmentation cost matrix on a single vector.
 *
 * Function re-implemented from Wolfgang Hubers tilingArray package
 * (https://github.com/Bioconductor-mirror/tilingArray/blob/master/R/costMatrix.R).
 * This is the version that runs on a single vector.
 *
 * @todo Profile this algorithm, to see where time is spent.
 * @todo Potentially replace calculate_cost_matrix_single by a faster version using
 *       boost UBLAS vectors and matrices.
 * 
 * @param data Single data vector.
 * @param max_k Maximum number of bins for each segment. max_k <= N.
 */
Matrix<double> calculate_cost_matrix(std::vector<double> const & data,
                                                int max_k)
{
    // See https://www.bioconductor.org/packages/devel/bioc/vignettes/tilingArray/inst/doc/costMatrix.pdf
    // for an explanation of the algebra

    double ZERO_THR = 1.1e-10;

    unsigned N = (unsigned)data.size();
    std::vector<double> cr(N);
    std::vector<double> cq(N);

    // cr = cumsum(data)
    std::partial_sum(data.begin(), data.end(), cr.begin());                             // O(N)

    // cq = cumsum(data^2)
    std::transform(data.begin(), data.end(), cq.begin(), [](double u){return u*u;});     // O(N)
    std::partial_sum(cq.begin(), cq.end(), cq.begin());                                 // O(N)

    // Cost matrix G
    Matrix<double> G(max_k, std::vector<double>(N, -1));                                  // O(N * maxk)

    // Initialization for segments of lengths 1...max_k at position 0.
    for (unsigned k = 0; k < max_k; ++k)                                                // O(maxk)
        G[k][0] = cq[k] - cr[k]*cr[k]/(double)(k+1);

    // Iterate per row (k)
    for (unsigned k = 1; k <= max_k && k <= N-1; ++k)                                   // O( maxk * ...
    {
        std::vector<double> cqk(N-k, -1);                                                //           ... N
        std::vector<double> crk(N-k, -1);                                                //           ... N

        for (unsigned i = 0; i < N-k; ++i) {                                            //           ... N
            cqk[i] = cq[i+k] - cq[i];
            crk[i] = cr[i+k] - cr[i];
        }

        // Iterate columns j = 1...N-k
        for (unsigned j = 1; j < N-k+1; ++j) {                                          //           ... N
            unsigned i = j - 1;
            assert(cqk[i]>=0);
            assert(crk[i]>=0);
            G[k-1][j] = cqk[i] - crk[i]*crk[i]/(double)k;

            // force small values to zero
            if (G[k-1][j] < ZERO_THR)
                G[k-1][j] = 0;
        }
    }

    return G;                                                                           // TOTAL: O(N * maxk) <= O(N^2)
}



/**
 * @fn calculate_cost_matrix(Matrix<double> const & data, int max_k, Matrix<double> & G)
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
bool calculate_cost_matrix(Matrix<double> const & data,                         // TOTAL O ( J * N * maxk)
                                    int max_k,
                                    Matrix<double> & G)
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


    G = Matrix<double>(max_k, std::vector<double>(N, 0));
    for (unsigned s = 0; s < J; ++s) {

        // Get cost of this one sample
        Matrix<double> cost_j = calculate_cost_matrix(data[s], max_k);

        // Add to total cost
        for (unsigned i=0; i<max_k; ++i)
            for (unsigned j=0; j<N; ++j)
                G[i][j] += cost_j[i][j];

    }

    // Finally, divide by J
    for (unsigned i=0; i<max_k; ++i)
        for (unsigned j=0; j<N; ++j)
            G[i][j] /= (double)J;

    return true;
}


struct Conf_segment {
    boost::filesystem::path f_in;
    boost::filesystem::path f_out;
    boost::filesystem::path f_cost_mat;
    float max_bp_per_Mb;
    float none_penalty;
    float merge_threshold;
    unsigned max_segment_length;
    bool remove_bad_cells;
};


/**
 * @fn main_segment(int argc, char** argv)
 * @ingroup segmentation
 * Main program for segmentation
 *
 * 1. Read Strand-seq counts from table, incl. 'class' column
 * 2. Remove cells which are all 'none' (opt-out).
 * 3. Detect stretches of removed bins (from 'None' in the counts table) and ...
 *    1. do nothing special,
 *    2. penalize these regions in the cost matrix, or
 *    3. remove these regions from the data before segmentation.
 * 4. Finally, run dynamic programming segmentation and report optimal segments
 *    for a different numbers of change points.
 */
int main_segment(int argc, char** argv) {

    Conf_segment conf;
    boost::program_options::options_description po_generic("Generic options");
    po_generic.add_options()
    ("help,?", "show help message")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("segments.txt"), "output file for counts")
    ;

    boost::program_options::options_description po_segmentation("Segmentation options");
    po_segmentation.add_options()
    ("max_bp,m", boost::program_options::value<float>(&conf.max_bp_per_Mb)->default_value(0.5), "maximum number of breakpoints per Mb")
    ("max_segment,M", boost::program_options::value<unsigned>(&conf.max_segment_length)->default_value(100000000), "maximum segment length")
    ("penalize-none", boost::program_options::value<float>(&conf.none_penalty)->implicit_value(100), "Penalize segments through removed bins (which are marked by 'None' in the counts table).")
    ("remove-none", "Remove segments through removed bins before segmentation. Mutually exclusive with --penalize-none.")
    ("do-not-remove-bad-cells", "Keep all cells (by default, cells which are marked 'None' in all bins get removed")
    ;

    boost::program_options::options_description po_hidden("Hidden options");
    po_hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&conf.f_in), "mosaicatcher count file")
    ("cost-matrix,c", boost::program_options::value<boost::filesystem::path>(&conf.f_cost_mat), "write cost matrix to file")
    ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", 1);


    boost::program_options::options_description cmdline_options;
    cmdline_options.add(po_generic).add(po_segmentation).add(po_hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(po_generic).add(po_segmentation);
    boost::program_options::variables_map vm;

    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);


    // Check arguments
    if (vm.count("help") || !vm.count("input-file")) {

        std::cout << std::endl;
        std::cout << "Mosaicatcher " << STRINGIFYMACRO(MOSAIC_VERSION_MAJOR);
        std::cout << "." << STRINGIFYMACRO(MOSAIC_VERSION_MINOR) << std::endl;
        std::cout << "> Find a segmentation across multiple cells." << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:   " << argv[0] << " [OPTIONS] counts.txt.gz" << std::endl << std::endl;
        std::cout << visible_options << std::endl;
        return vm.count("help") ? 0 : 1;
    }
    conf.remove_bad_cells = true;
    if (vm.count("do-not-remove-bad-cells"))
        conf.remove_bad_cells = false;
    if (vm.count("penalize-none") && vm.count("remove-none")) {
        std::cerr << "[Error] --penalize-none and --remove-none are mutually exclusive." << std::endl;
        return 1;
    }



    /////////////////////////////////////////////////////////// global variables
    /* counts/cells */
    std::vector<std::pair<std::string,std::string>> sample_cell_names;
    std::vector<TGenomeCounts> counts;

    /* bins */
    std::vector<std::string> chromosomes;
    std::vector<Interval> bins;
    std::vector<int32_t>  chrom_map;
    std::vector<unsigned> good_bins;
    std::vector<int32_t>  good_map;
    /////////////////////////////////////////////////////////// global variables


    std::chrono::steady_clock::time_point t1, t2, t3, t4;


    // 1. Reading count file
    t1 = std::chrono::steady_clock::now();
    std::cout << "[Info] Reading file " << conf.f_in.string() << " (this is a bit slow)..." << std::endl;
    if (!io::read_counts_gzip(conf.f_in.string(),
                              counts,
                              chromosomes,
                              sample_cell_names,
                              bins))
        return 1;

    t2 = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "[Info] " << counts.size() << " cells found." << std::endl;
    std::cout << "[Time] Reading took " << time.count() << " seconds." << std::endl;

    // If more than one sample - give a warning
    std::set<std::string> samples;
    for (auto s_c : sample_cell_names) samples.insert(s_c.first);
    if (samples.size()>1)
        std::cerr << "[Warning] Found > 1 sample, but Samples will be ignored for now. All cells will be treated as being from the first sample" << std::endl;


    // Create chromosome map
    chrom_map = std::vector<int32_t>(chromosomes.size());
    if (!make_chrom_map(bins, chrom_map)) {
        std::cerr << "[Error] Something is wrong with the intervals" << std::endl;
        return 2;
    }
    chrom_map.push_back((int32_t)bins.size());
    std::cout << "[Info] " << bins.size() << " intervals across " << chromosomes.size() << " chromosomes." << std::endl;


    // 2. remove bad cells (i.e. with all None).
    std::vector<unsigned> good_cells = io::get_good_cells(counts);
    if (conf.remove_bad_cells) {
        unsigned prev = counts.size();

        for (unsigned i = 0; i < good_cells.size(); ++i) {
            counts[i] = counts[good_cells[i]];
        }
        counts.resize(good_cells.size()); // should
        std::cout << "[Info] Removed " << prev - counts.size() << " bad cells, leaving " << counts.size() << " good ones" << std::endl;
    }



    // Determine window size & Mean count per sample
    unsigned window_size;
    std::vector<float> mean_per_sample(counts.size());
    {
        TMeanVarAccumulator<float> mean_acc;
        for (auto bin : bins)
            mean_acc((float)(bin.end - bin.start));
        window_size = (unsigned) boost::accumulators::mean(mean_acc);

        for (unsigned i=0; i<counts.size(); ++i) {
            TMeanVarAccumulator<float> mean_acc;
            for (auto j = 0; j < bins.size(); ++j)
                mean_acc((counts[i][j]).watson_count + (counts[i][j]).crick_count);
            mean_per_sample[i] = boost::accumulators::mean(mean_acc);
        }
    }





    // Determine `good_bins` from stretches of 'None' in the data.
    // This is more flexible than using the count data, because it allows
    // the user to alter the labels of certain regions!
    std::vector<std::pair<unsigned,unsigned>> stretches;
    if (vm.count("penalize-none") || vm.count("remove-none"))
    {
        // Finding 'none' stretches
        for (int32_t chrom = 0; chrom < chromosomes.size(); ++chrom)
        {
            unsigned pos = chrom_map[chrom];
            while (pos < chrom_map[chrom+1]) {
                if (counts[0][pos].label != "None") {
                    good_bins.push_back(pos);
                    ++pos;
                    continue;
                }
                // iterate through consecutive stretch of None
                unsigned start = pos;
                while (counts[0][pos].label == "None" && pos < chrom_map[chrom+1])
                    ++pos;
                stretches.push_back(std::make_pair(start, pos-1));
                std::cout << "    > Discovered 'None' stretch in bins [" << start << "," << pos -1 << "]" << std::endl;
            }

            // build chrom_map for good bins
            good_map = std::vector<int32_t>(chrom_map.size() - 1, -1);
            pos = 0;
            for (int32_t chr = 0; chr < static_cast<int32_t>(good_map.size()); ++chr) {
                while (pos < good_bins.size() && bins[good_bins[pos]].chr < chr)
                    ++pos;
                // now goodit is either at first occurence of chr, or at the end.
                if (pos >= good_bins.size()) good_map[chr] = (int32_t)good_bins.size();
                else good_map[chr] = pos;
            }
            // add last element for easy calculation of number of bins
            good_map.push_back((int32_t)good_bins.size());
        }
    } else {
        good_bins.resize(bins.size());
        std::iota(good_bins.begin(), good_bins.end(), 0); // 0,1,2,...
        good_map = chrom_map;
    }







    // prepare OUTPUT file
    std::ofstream out(conf.f_out.string());
    out << "sample" << "\t" << "cells" << "\t" << "chrom" << "\t" << "bins";
    out << "\t" << "maxcp" << "\t" << "none" << "\t" << "action" << "\t" << "k";
    out << "\t" << "sse" << "\t" << "bps" << "\t" << "start" << "\t" << "end";
    out << std::endl;


    // Segmentation chromosome per chromosome
    // This loop can later be parallelized
    for (int32_t chrom=0; chrom < chromosomes.size(); ++chrom)
    {

        unsigned N = chrom_map[chrom+1] - chrom_map[chrom];
        if (N < 1) continue;

        // parameters:
        unsigned chrom_size = bins[chrom_map[chrom+1]-1].end - bins[chrom_map[chrom]].start;
        // max_k = longest allowed segment = 20Mb or max chrom_size
        unsigned max_k  = std::min(std::max(static_cast<unsigned>(10),
                                            static_cast<unsigned>(conf.max_segment_length/window_size)),
                                   static_cast<unsigned>(N));
        // max_cp = max. number of change points = max_bp_per_Mb * Mb, but at least 10
        unsigned max_cp = std::min(std::max(static_cast<unsigned>(10),
                                            static_cast<unsigned>(ceil((float)chrom_size/1e6 * conf.max_bp_per_Mb))),
                                   N-1);


        // Put all cells into data:
        t1 = std::chrono::steady_clock::now();
        Matrix<double> data;
        // temporarily write counts into this vector before adding to data
        std::vector<double> tmp(N);
        for (unsigned i = 0; i < counts.size(); ++i)
        {
            double sample_mean = mean_per_sample[i];

            std::transform(counts[i].begin() + chrom_map[chrom],
                           counts[i].begin() + chrom_map[chrom] + N,
                           tmp.begin(),
                           [sample_mean](Counter const & c){return (double) c.watson_count / sample_mean;});
            data.push_back(tmp);
            std::transform(counts[i].begin() + chrom_map[chrom],
                           counts[i].begin() + chrom_map[chrom] + N,
                           tmp.begin(),
                           [sample_mean](Counter const & c){return (double) c.crick_count / sample_mean;});
            data.push_back(tmp);
        }
        tmp.clear();
        t2 = std::chrono::steady_clock::now();
        std::cout << "    Preparing `data` took " << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << " seconds." << std::endl;



        // Mode 2
        // Remove 'none' bins from data before running segmentation.
        // Later I will have to add these gaps again!
        int num_none = 0;
        if (vm.count("remove-none"))
        {
            std::cout << "    Removing None stretches" << std::endl;
            num_none = (chrom_map[chrom+1] - chrom_map[chrom]) - (good_map[chrom+1] - good_map[chrom]);
            if (num_none == 0) {
                std::cout << "    No bad bins found." << std::endl;
            } else {
                // Recalculate N, max_cp, and max_k
                N = data[0].size() - num_none;
                max_k  = std::min(max_k, N);
                max_cp = std::min(max_cp, N);
                // Remove bins (columns) from data matrix
                for (unsigned i = 0; i < data.size(); ++i) {
                    for (unsigned gbin = good_map[chrom], lbin = 0; gbin < good_map[chrom+1]; ++gbin, ++lbin) {
                        data[i][lbin] = data[i][good_bins[gbin]];
                    }
                    data[i].resize(N);
                }
                std::cout << "    Removed " << num_none << " bins" << std::endl;
            }
        }


        // New Cost matrix
        t1 = std::chrono::steady_clock::now();
        Matrix<double> new_cost;
        calculate_cost_matrix(data, max_k, new_cost);

        // also, check that values make sense
        for (unsigned k=0; k<max_k; ++k)
            for (unsigned j=0; j<N-k; ++j)
                assert(new_cost[k][j] >= 0);
        for (unsigned k=0; k<max_k; ++k)
            for (unsigned j=N-k; j<N; ++j)
                assert(new_cost[k][j] <= 0);




        // Mode 1 (none bins)
        // Force segmentation to use "None" intervals.
        // This is done by setting the cost for such intervals to 0 and adding
        // a penalty to intervals violating the boarders.

        if (vm.count("penalize-none"))
        {
            std::cout << "    Penalty for None = " << conf.none_penalty << std::endl;
            for (auto stretch : stretches)
            {
                // Penalize all stretches that violate these boarders.
                // Todo: note that this penalty becomes weaker the longer the stretches are!
                for (unsigned k = 0; k < max_k; ++k)
                    for (unsigned j = (k>stretch.first ? 0 : stretch.first - k); j < stretch.second && j < N-k; ++j)
                        new_cost[k][j] += conf.none_penalty;

                // Finally, favor the usage of exactly the consecutive None segment
                new_cost[stretch.second - stretch.first  - 1][stretch.first ] = 0;
            }
        } // vm.count("penalize-none")


        t2 = std::chrono::steady_clock::now();
        std::cout << "    Cost matrix (" << max_k << " x " << N << ") took " << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << " seconds." << std::endl;



        // print cost matrix for first chromosome
        // todo: remove later
        if (vm.count("cost-mat") && chrom == 0) {
            std::ofstream out(conf.f_cost_mat.string());
            if (out.is_open()) {
                std::cout << "[Write] cost matrix for 1. chromosome " << conf.f_cost_mat.string() << std::endl;
                out << new_cost;
            } else {
                std::cerr << "[Warning] Cannot write to " << conf.f_cost_mat.string() << std::endl;
            }
        }



        // Find optimal segmentation
        Matrix<int> breakpoints;
        std::vector<double> sse;
        if (!optimal_segment_dp(new_cost, max_cp, breakpoints, sse))
            std::cerr << "[ERROR] Segmentation failed" << std::endl;



        // Output of breakpoints;
        if (out.is_open()) {
            std::cout << "    [Write] " << chromosomes[chrom] << " to file: " << conf.f_out.string() << std::endl;
            for (unsigned cp = 0; cp < max_cp; ++cp) {
                for (unsigned k = 0; k <= cp; ++k) {
                    out << *samples.begin() <<  "\t";
                    out << counts.size() << "\t";
                    out << chromosomes[chrom] << "\t";
                    out << chrom_map[chrom+1] - chrom_map[chrom] << "\t";
                    out << max_cp << "\t";
                    out << num_none << "\t";
                    if (vm.count("penalize-none"))
                        out << "penalize-none" << "\t";
                    else if (vm.count("remove-none"))
                        out << "remove-none" << "\t";
                    else
                        out << "utilize-none" << "\t";
                    out << cp+1 << "\t";
                    out << sse[cp] << "\t";
                    unsigned from = k==0 ? 0 : breakpoints[cp][k-1];
                    unsigned to   = breakpoints[cp][k]-1;
                    if (vm.count("remove-none")) {
                        from = good_bins[from];
                        to   = good_bins[to];
                    }
                    out << to << "\t";
                    out << bins[from].start << "\t";
                    out << bins[to].end << std::endl;
                }
            }
        } else {
            std::cerr << "[Warning] Cannot write to " << conf.f_out.string() << std::endl;
        }

    } // for chrom

    return 0;
} // main



