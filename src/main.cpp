/*
Copyright (C) 2016 Sascha Meiers
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Sascha Meiers (meiers@embl.de)
============================================================================
*/

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <htslib/sam.h>

#include "intervals.hpp"
#include "counts.hpp"
#include "distribution.hpp"
#include "hmm.hpp"


static unsigned MIN_MEDIAN_PER_SAMPLE = 20;
static double   MIN_MEDIAN_PER_BIN    = 0.1;


struct Conf {
    std::vector<boost::filesystem::path> f_in;
    boost::filesystem::path f_out;
    boost::filesystem::path f_bins;
    boost::filesystem::path f_excl;
    int minMapQual;
    unsigned int window;
    std::string mode;
};

struct CellInfo {
    unsigned median_bin_count;
    std::string sample_name;
};

std::vector<unsigned> median_by_sample(std::vector<TGenomeCounts> & counts)
{
    std::vector<unsigned> median_by_sample(counts.size());
    for(unsigned i = 0; i < counts.size(); ++i) {
        TMedianAccumulator<unsigned int> med_acc;
        for (Counter const & count_bin : counts[i])
            med_acc(count_bin.watson_count + count_bin.crick_count);
        median_by_sample[i] = boost::accumulators::median(med_acc);
    }
    return median_by_sample;
}



std::vector<double> median_per_bin(std::vector<TGenomeCounts> & counts)
{
    std::vector<double> median_per_bin(counts[0].size());
    for (unsigned bin = 0; bin < counts[0].size(); ++bin)
    {
        TMedianAccumulator<double> med_acc;
        for(unsigned i = 0; i < counts.size(); ++i) {
            Counter const & cc = counts[i][bin];
            med_acc(cc.watson_norm + cc.crick_norm);
        }
        median_per_bin[bin] = boost::accumulators::median(med_acc);
    }
    return median_per_bin;
}



void run_gaussian_HMM(std::vector<TGenomeCounts> & counts, std::vector<int32_t> const & chrom_map)
{
    /*
     hmm::MultiVariateGaussianHMM hmm(12, 2);
     std::vector<double> transition_probs = {\
     981,  5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
     5,  981, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
     5,  5, 981, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
     5,  5, 5, 977, 1, 1, 1, 1, 1, 1, 1, 1, \
     5,  5, 5, 1, 977, 1, 1, 1, 1, 1, 1, 1, \
     5,  5, 5, 1, 1, 977, 1, 1, 1, 1, 1, 1, \
     5,  5, 5, 1, 1, 1, 977, 1, 1, 1, 1, 1, \
     5,  5, 5, 1, 1, 1, 1, 977, 1, 1, 1, 1, \
     5,  5, 5, 1, 1, 1, 1, 1, 977, 1, 1, 1, \
     5,  5, 5, 1, 1, 1, 1, 1, 1, 977, 1, 1, \
     5,  5, 5, 1, 1, 1, 1, 1, 1, 1, 977, 1, \
     5,  5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 977 };
     for (unsigned i=0; i< transition_probs.size(); ++i)
     transition_probs[i] /= 1000;
     hmm.set_transitions(transition_probs);
     hmm.setInitials({0.25,0.41,0.25, 0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01 });
     hmm.set_emissions( {hmm::Gaussian("WW",   2, {0,2},   {0.1,0.4}),
     hmm::Gaussian("WC",   2, {1,1},   {0.4,0.4}),
     hmm::Gaussian("CC",   2, {2,0},   {0.4,0.1}),
     hmm::Gaussian("CCC",  2, {3,0},   {0.4,0.1}),
     hmm::Gaussian("WCC",  2, {2,1},   {0.4,0.1}),
     hmm::Gaussian("WWC",  2, {1,2},   {0.4,0.1}),
     hmm::Gaussian("WWW",  2, {0,3},   {0.4,0.1}),
     hmm::Gaussian("CCCC", 2, {4,0},   {0.4,0.1}),
     hmm::Gaussian("CCCW", 2, {3,1},   {0.4,0.1}),
     hmm::Gaussian("CCWW", 2, {2,2},   {0.4,0.1}),
     hmm::Gaussian("CWWW", 2, {1,3},   {0.4,0.1}),
     hmm::Gaussian("WWWW", 2, {0,4},   {0.4,0.1})    } );
     */



    // Set up an HMM:
    hmm::HMM<double, hmm::MultiVariate<hmm::Gaussian> > hmm({"WW", "WC", "CC"});
    hmm.set_transitions({\
        .998,  .001, .001, \
        .001, .998, .001, \
        .001, .001, .998 });
    hmm.set_emissions( {\
        hmm::MultiVariate<hmm::Gaussian>({hmm::Gaussian(0,0.5), hmm::Gaussian(2,0.5)}),
        hmm::MultiVariate<hmm::Gaussian>({hmm::Gaussian(1,0.5), hmm::Gaussian(1,0.5)}),
        hmm::MultiVariate<hmm::Gaussian>({hmm::Gaussian(2,0.7), hmm::Gaussian(0,0.5)})
    });


    for (unsigned i = 0; i < counts.size(); ++i) {

        // Order: crick, watson, crick, watson, ...
        std::vector<double> seq;
        unsigned chr_idx = 0;
        for (unsigned bin = 0; bin < counts.size(); ++bin) {
            seq.push_back(counts[i][bin].crick_norm);
            seq.push_back(counts[i][bin].watson_norm);

            // chromosome finished:
            if (bin == chrom_map[chr_idx]) {
                if (seq.size()>0) {

                    // run HMM
                    hmm.viterbi(seq);
                    std::vector<std::string> path = hmm.get_path_labels();
                    assert(path.size() == seq.size()/2);

                    // write classification into Counter
                    unsigned bin_in_path = 0;
                    for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
                        Counter & cc = counts[i][bin];
                        cc.set_label(path[bin_in_path++]);
                    }
                }
                seq.clear();
            }
        }
    }
}


void run_bn_HMM(std::vector<TGenomeCounts> & counts, std::vector<int32_t> const & chrom_map, std::vector<unsigned> const & sample_median)
{
    assert(counts.size() == sample_median.size());

    // Set up an HMM:
    hmm::HMM<unsigned, hmm::MultiVariate<hmm::NegativeBinomial> > hmm({"WW", "WC", "CC"});
    hmm.set_transitions({\
        .998,  .001, .001, \
        .001, .998, .001, \
        .001, .001, .998 });

    for (unsigned i=0; i<counts.size(); ++i) {

        // set n and p parameters according to a mean of "sample_median"
        // n = 1
        double p = 0.2;
        double med = (double)sample_median[i]/2;
        double n = (double)med * p / (1-p);
        hmm.set_emissions( {\
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,2*n), hmm::NegativeBinomial(p,  1)}),
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  n), hmm::NegativeBinomial(p,  n)}),
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  1), hmm::NegativeBinomial(p,2*n)})
        });


        { // print
                            std::cout << std::endl;
                            std::cout << "----------------------------" << std::endl;
                            std::cout << "Sample " << i << "\t" << "mean = " << med << "\t" << "n = " << n << "\t" << "p = " << p << std::endl;

                            for (auto dist : hmm.distributions) {
                                std::cout << dist;
                            }
        } // end print


        for (unsigned i = 0; i < counts.size(); ++i) {

            // Order: crick, watson, crick, watson, ...
            std::vector<unsigned> seq;
            unsigned chr_idx = 0;
            for (unsigned bin = 0; bin < counts.size(); ++bin) {
                seq.push_back(counts[i][bin].crick_count);
                seq.push_back(counts[i][bin].watson_count);

                // chromosome finished:
                if (bin == chrom_map[chr_idx]) {
                    if (seq.size()>0) {

                        // run HMM
                        hmm.viterbi(seq);
                        std::vector<std::string> path = hmm.get_path_labels();
                        assert(path.size() == seq.size()/2);

                        // write classification into Counter
                        unsigned bin_in_path = 0;
                        for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
                            Counter & cc = counts[i][bin];
                            cc.set_label(path[bin_in_path++]);
                        }
                    }
                    seq.clear();
                }
            }
        }
    }
}


int main(int argc, char **argv)
{

    // Command line options
    Conf conf;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
    ("help,?", "show help message")
    ("mapq,q", boost::program_options::value<int>(&conf.minMapQual)->default_value(10), "min mapping quality")
    ("window,w", boost::program_options::value<unsigned int>(&conf.window)->default_value(1000000), "window size of fixed windows")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("out.txt"), "output file for counts")
    ("bins,b", boost::program_options::value<boost::filesystem::path>(&conf.f_bins), "variable bin file (BED format, mutually exclusive to -w)")
    ("exclude,x", boost::program_options::value<boost::filesystem::path>(&conf.f_excl), "Exclude chromosomes (mutually exclusive to -b)")
    ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
    ("input-file", boost::program_options::value<std::vector<boost::filesystem::path> >(&conf.f_in), "input bam file(s)")
    ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;

    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if (!vm["window"].defaulted() && vm.count("bins")) {
        std::cerr << "Error: -w and -b cannot be specified together" << std::endl << std::endl;
        goto print_usage_and_exit;
    }
    if (vm.count("bins") && vm.count("exclude")) {
        std::cerr << "Error: Exclude chromosomes (-x) have no effect when -b is specifiet. Stop" << std::endl << std::endl;
        goto print_usage_and_exit;
    }
    if (vm.count("help") || !vm.count("input-file"))
    {
    print_usage_and_exit:
        std::cout << "Usage: " << argv[0] << " [OPTIONS] <strand.seq1.bam> <strand.seq2.bam> ... <strand.seqN.bam>" << std::endl;
        std::cout << visible_options << std::endl;
        std::cout << std::endl;
        std::cout << "Notes:" << std::endl;
        std::cout << "  * Reads are counted by start position" << std::endl;
        std::cout << "  * One cell per BAM file, inclusing SM tag in header" << std::endl;
        std::cout << "  * For paired-end data, only read 1 is counted" << std::endl;
        return 1;
    }


    // global variables
    bam_hdr_t* hdr = NULL;
    std::vector<CellInfo>      cells(conf.f_in.size());
    std::vector<Interval>      bins;
    std::vector<int32_t>       chrom_map;
    std::vector<TGenomeCounts> counts(conf.f_in.size());



    // Read sample names from headers.
    // Keep one header throughout the program.
    std::cout << "Reading SAM headers" << std::endl;
    for(int i = 0; i < conf.f_in.size(); ++i)
    {
        samFile* samfile = sam_open(conf.f_in[i].string().c_str(), "r");
        if (samfile == NULL) {
            std::cerr << "Fail to open file " << conf.f_in[0].string() << std::endl;
            return 1;
        }
        hdr = sam_hdr_read(samfile);
        if (!get_SM_tag(hdr->text, cells[i].sample_name)) {
            std::cerr << "Each BAM file has to have exactly one SM tag." << std::endl << std::endl;
            goto print_usage_and_exit;
        }
    }



    // Bin the genome
    chrom_map = std::vector<int32_t>(hdr->n_targets, -1);
    if (vm.count("bins"))
    {
        if (!read_dynamic_bins(bins, chrom_map, conf.f_bins.string().c_str(), hdr))
            return 1;
    }
    else
    {
        std::vector<Interval> exclude;
        if (vm.count("exclude")) {
            read_exclude_file(conf.f_excl.string(), hdr, exclude);
            sort(exclude.begin(), exclude.end(), interval_comp);
        }
        std::cout << "Exclude " << exclude.size() << " regions" << std::endl;
        create_fixed_bins(bins, chrom_map, conf.window, exclude, hdr);
    }
    // add last element for easy calculation of number of bins
    chrom_map.push_back((int32_t)bins.size());



    // Count in bins
    for(int i = 0; i < counts.size(); ++i)
    {
        if (!count_sorted_reads(conf.f_in[i].string(), bins, chrom_map, hdr, conf.minMapQual, counts[i])) {
            std::cerr << "Ignore sample " << conf.f_in[i].string() << std::endl;
            counts.erase(counts.begin()+i);
            --i;
        }

        // todo: kick out samples if median is too low
    }

    // Mode
    std::cout << "Writing " << conf.f_out.string() << std::endl;
    std::ofstream out(conf.f_out.string());
    out << "chrom\tstart\tend\tsample\tw\tc" << std::endl;
    for(unsigned i = 0; i < counts.size(); ++i) {
        for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
            out << hdr->target_name[bins[bin].chr];
            out << "\t" << bins[bin].start << "\t" << bins[bin].end;
            out << "\t" << conf.f_in[i].filename().string();
            out << "\t" << counts[i][bin].watson_count;
            out << "\t" << counts[i][bin].crick_count;
            out << std::endl;
        }
    }
    out.close();

    // todo: string comparisons as labels --> replace by faster type (e.g. enum)





    std::cout << "Caution: Classification of strandedness is still very preliminary" << std::endl;



    // 1. Normalize counts by sample:
    //    To do: How to handle samples that are kicked out?
    std::vector<unsigned> sample_median = median_by_sample(counts);
    for(unsigned i = 0; i < counts.size(); ++i) {

        if (sample_median[i] < MIN_MEDIAN_PER_SAMPLE) {
            std::cout << "    Todo(!): Ignoring sample " << i << std::endl;
            continue;
        }
        for(Counter & cc : counts[i]) {
            cc.watson_norm = cc.watson_count / (double)sample_median[i] * 2;
            cc.crick_norm  = cc.crick_count  / (double)sample_median[i] * 2;
        }
    }


    // 3. Run HMM
    bool gauss=false;
    if (gauss)
        run_gaussian_HMM(counts, chrom_map);
    else
        run_bn_HMM(counts, chrom_map, sample_median);



    // print normalized counts
    if (conf.mode == "classify")
    {
        std::cout << "Writing file " << conf.f_out.string() << std::endl;
        std::ofstream out(conf.f_out.string());
        out << "chrom\tstart\tend\tsample\tw\tc\twn\tcn\tclass" << std::endl;
        for(unsigned i = 0; i < counts.size(); ++i) {
            unsigned chrom = 0;
            for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
                Counter & cc = counts[i][bin];
                out << hdr->target_name[bins[bin].chr];
                out << "\t" << bins[bin].start << "\t" << bins[bin].end;
                out << "\t" << conf.f_in[i].filename().string();
                out << "\t" << cc.watson_count;
                out << "\t" << cc.crick_count;
                out << "\t" << cc.watson_norm;
                out << "\t" << cc.crick_norm;
                out << "\t" << cc.get_label();
                out << std::endl;
            }
        }
        out.close();
        return(0);
    } 


    // nothing here

    return(0);
}
