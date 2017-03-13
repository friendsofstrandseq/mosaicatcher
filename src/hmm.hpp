//
//  hmm.hpp
//  hmm
//
//  Created by Sascha Meiers on 24/02/2017.
//  Copyright © 2017 Sascha Meiers. All rights reserved.
//

#ifndef hmm_hpp
#define hmm_hpp

#include <iostream>
#include <vector>
#include <cassert>
#include "utils.hpp"
#include <boost/math/distributions/normal.hpp>


// todo: remove if std::setprecision no longer needed
#include <iomanip>

namespace hmm {

    static const double INPUT_PRECISION = 0.0001;

    /**
     *  Gaussian
     *  --------
     *  stores parameters for a multivariate normal distribution
     *  Unfortunately boost::math::normal does not support a covariance matrix,
     *  so components must be independent (covariance = 0)
     */
    struct Gaussian {
        std::string label;
        uint8_t dim;
        std::vector<double> mean;
        std::vector<double> var;

        Gaussian(std::string const & label,
                 uint8_t dim,
                 std::vector<double> const & mean,
                 std::vector<double> const & var) : label(label), dim(dim) {
            assert(mean.size() == dim);
            assert(var.size() == dim);
            this->mean = mean;  // copy
            this->var  = var;
        }
    };



    // Note: Do not copy this class (no copy constructor defined)
    class MultiVariateGaussianHMM {
    public:
        // Initialization:
        MultiVariateGaussianHMM(uint8_t N, uint8_t d);
        ~MultiVariateGaussianHMM();
        void set_transitions(std::vector<double> const &);
        void setInitials(std::vector<double> const &);
        void set_emissions(std::vector<Gaussian> const &);

        // Access to variables/results
        std::vector<uint8_t> const & get_path() {return path; }
        std::vector<std::string> get_path_labels();

        // Algorithms
        double viterbi(std::vector<double> const &); // writes "path"
        std::vector<std::vector<double> > forward(std::vector<double> const &);
        //std::vector<std::vector<double> > backward(std::vector<double> const &);

    private:
        std::vector<double> calc_log_emissions(std::vector<double>::const_iterator iter);
        uint8_t N;           // Number of states
        uint8_t d;           // how many emissions per time step
        double* pi;          // initial transition probabilites
        double** log_trans;  // transition probabilites (log scale)
        std::vector<Gaussian> emission_params;   // emission parameters - one "Gaussian" per state
        std::vector<uint8_t> path;      // last path calculated;
    };

    MultiVariateGaussianHMM::MultiVariateGaussianHMM(uint8_t N, uint8_t d) :
    N(N), d(d)
    {
        pi = new double[N];
        log_trans = new double*[N];
        for (uint8_t i=0; i<N; ++i)
            log_trans[i] = new double[N];
    }

    MultiVariateGaussianHMM::~MultiVariateGaussianHMM()
    {
        delete[] pi;
        for (uint8_t i=0; i<N; ++i)
            delete[] log_trans[i];
        delete[] log_trans;
    }

    void MultiVariateGaussianHMM::set_transitions(std::vector<double> const & t)
    {
        assert(t.size() == N*N);
        for (uint8_t i=0; i<N; ++i) {
            double sum = 0;
            for (uint8_t j=0; j<N; ++j) {
                sum += t[N*i+j];
                log_trans[i][j] = log(t[N*i+j]);
            }
            assert(std::abs(1 - sum) < INPUT_PRECISION);
        }
    }

    void MultiVariateGaussianHMM::setInitials(std::vector<double> const & ini) {
        assert(ini.size() == N);
        assert(sum(ini) == 1.0);
        for (uint8_t i=0; i<N; ++i)
            pi[i] = ini[i];
    }

    void MultiVariateGaussianHMM::set_emissions(std::vector<Gaussian> const & e)
    {
        assert(e.size() == N);
        for (Gaussian const & x : e) {
            assert(x.dim == d);
        }
        emission_params = e; // copy!
    }

    std::vector<std::string> MultiVariateGaussianHMM::get_path_labels() {
        std::vector<std::string> labels;
        for (unsigned i=0; i<path.size(); i++) {
            labels.push_back(emission_params[path[i]].label);
        }
        return labels;
    }


    /**
     *  Given an observation (which consists of multiple doubles for a 
     *  multivariate model), calculate the PDFs for all states.
     *  Result should be move-constructed efficiently
     *  TODO: how can I whether it is really moved and not copied?
     */
    std::vector<double> MultiVariateGaussianHMM::calc_log_emissions(std::vector<double>::const_iterator iter)
    {
        std::vector<double>::const_iterator start = iter;
        std::vector<double> log_probs(N);

        for (uint8_t state = 0; state < N; ++state)
        {
            Gaussian const& par = emission_params[state];
            double log_prob = 0;
            // for each component, sum up their log probabilites.
            // multivariate emissions are encoded succeedingly in one long vector
            for (uint8_t j=0; j<d; ++j, ++iter) // iter could run over sequence, careful!
                log_prob += log(pdf(boost::math::normal(par.mean[j], par.var[j]), *iter));
            log_probs[state] = log_prob;
            iter = start;
        }
        return log_probs;
    }


    /**
     *  Viterbi algorithm
     *  -----------------
     *  Sequence must be flattened for multivariate emissions, i.e. a 2-dim
     *  obervation across T=4 time points, e.g. [(0,1), (2,4), (1,3), (4,5)],
     *  would be written [0,1,2,4,1,3,4,5].
     */
    double MultiVariateGaussianHMM::viterbi(std::vector<double> const & seq)
    {
        // Nothing to do for empty sequences
        if (seq.size() < d) {
            path = std::vector<uint8_t>();
            return 0;
        }

        // length of the sequence
        unsigned T = (unsigned)seq.size()/d;

        // allocate a matrix for probabilities:
        // matrix[T=columns][N=rows]
        std::vector<std::vector<double> > matrix(T, std::vector<double>(N, 1));

        // allocate a matrix for the path through the matrix
        std::vector<std::vector<uint8_t> > trace(T, std::vector<uint8_t>(N, N));

        // Step 1: Initiation
        // Prepare emissions for first observation:
        std::vector<double>::const_iterator it = seq.begin();
        std::vector<double> log_emission = calc_log_emissions(it);
        it += d; // First observation was already processed

        for (uint8_t state = 0; state < N; ++state)
            matrix[0][state] = log(pi[state]) + log_emission[state];


        // Step 2: Iteration
        for (unsigned i = 1; i < T; ++i, it+=d)
        {
            // get emission for observation(s) i
            log_emission = calc_log_emissions(it);


            for (uint8_t state = 0; state < N; ++state)
            {
                double max_log = -10000; // todo: must actually be smaller than all log(p)
                uint8_t max_path = N;

                for (uint8_t prev_state = 0; prev_state < N; ++prev_state) {
                    double log_prob = log_trans[prev_state][state] + matrix[i-1][prev_state];
                    if (log_prob > max_log) {
                        max_log = log_prob;
                        max_path = prev_state;
                    }
                }
                matrix[i][state] = max_log + log_emission[state];
                trace[i][state] = max_path;
            }
        }


        // Get maximum in last column
        double max_log = -10000;
        uint8_t max_path = N;
        for (uint8_t state = 0; state < N; ++state) {
            if (matrix[T-1][state] > max_log) {
                max_log = matrix[T-1][state];
                max_path = state;
            }
        }

        //std::cout << "max path " << (int)max_path << std::endl;

        // Step 3: Traceback
        // overwrite the member "path"
        path = std::vector<uint8_t>(T, N);
        for (unsigned i=T; i>0; ) {
            --i;
            path[i] = max_path;
            max_path = trace[i][max_path];
        }

        bool print_matrix = false;
        if (print_matrix)
        { // print

            it = seq.begin();
            log_emission = calc_log_emissions(it);
            it += d; // First observation was already processed


            // header
            std::cout << std::endl;
            std::cout << std::setw(6) << "bin    matrix                          path                            emission                        label" << std::endl;
            std::cout << std::setw(6) << "number";
            for (uint8_t state = 0; state < N; ++state)
                std::cout << std::setw(8) << (int) state;
            for (uint8_t state = 0; state < N; ++state)
                std::cout << std::setw(8) << (int) state;
            for (uint8_t state = 0; state < N; ++state)
                std::cout << std::setw(8) << (int) state;
            std::cout << std::endl;
            std::cout << "------ ------------------------------- ------------------------------- ------------------------------- --------" << std::endl;

            // matrix and path
            for (unsigned i = 0; i < T; ++i)
            {

                log_emission = calc_log_emissions(it);
                it += d;

                std::cout << std::setw(6) << (int)i;

                // first cols: matrix
                for (uint8_t state = 0; state < N; ++state)
                    std::cout << std::setw(8) << std::setw(8) << std::setprecision(1) << matrix[i][state];

                // second set of cols: path
                for (uint8_t state = 0; state < N; ++state)
                    std::cout << std::setw(8) << (int)trace[i][state];

                for (uint8_t state = 0; state < N; ++state)
                    std::cout << std::setw(8) << std::setprecision(2) << log_emission[state];

                std::cout << std::setw(8) << std::setprecision(2) << emission_params[path[i]].label;

                std::cout << std::endl;

            }

        } // end of print

        return max_log;
    }


    /**
     *  Forward algorithm
     */
    std::vector<std::vector<double> > MultiVariateGaussianHMM::forward(std::vector<double> const & seq)
    {
        // length of the sequence
        unsigned T = (unsigned)seq.size()/d;

        // allocate a matrix for probabilities:
        std::vector<std::vector<double> > matrix(T, std::vector<double>(N, 1));

        // Step 1: Initiation
        // Prepare emissions for first observation:
        std::vector<double>::const_iterator it = seq.begin();
        std::vector<double> log_emission = calc_log_emissions(it);
        it += d; // First observation was already processed

        for (uint8_t state = 0; state < N; ++state)
            matrix[0][state] = log(pi[state]) + log_emission[state];

/*
        // Step 2: Iteration
        for (unsigned i = 1; i < T; ++i, it+=d)
        {
            // get emission for observation(s) i
            log_emission = calc_log_emissions(it);

            for (uint8_t state = 0; state < N; ++state)
            {
                double max_log = -10000; // todo: must actually be smaller than all log(p)
                uint8_t max_path = N;

                for (uint8_t prev_state = 0; prev_state < N; ++prev_state) {
                    double log_prob = log_trans[prev_state][state] + matrix[i-1][prev_state];
                    if (log_prob > max_log) {
                        max_log = log_prob;
                        max_path = prev_state;
                    }
                }
                matrix[i][state] = max_log + log_emission[state];
            }
        }
*/

        return matrix;
    }


}
#endif /* hmm_hpp */
