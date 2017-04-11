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


struct Conf {
    std::vector<boost::filesystem::path> f_in;
    boost::filesystem::path f_out;
    boost::filesystem::path f_bins;
    boost::filesystem::path f_excl;
    boost::filesystem::path f_info;
    int minMapQual;
    unsigned int window;
    std::string mode;
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





int main(int argc, char **argv)
{

    // Command line options
    Conf conf;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
    ("help,?", "show help message")
    ("mode,m", boost::program_options::value<std::string>(&conf.mode)->default_value("hmm"), "mode: count | hmm")
    ("mapq,q", boost::program_options::value<int>(&conf.minMapQual)->default_value(10), "min mapping quality")
    ("window,w", boost::program_options::value<unsigned int>(&conf.window)->default_value(1000000), "window size of fixed windows")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("out.txt"), "output file for counts")
    ("bins,b", boost::program_options::value<boost::filesystem::path>(&conf.f_bins), "variable bin file (BED format, mutually exclusive to -w)")
    ("exclude,x", boost::program_options::value<boost::filesystem::path>(&conf.f_excl), "Exclude chromosomes (mutually exclusive to -b)")
    ("info,i", boost::program_options::value<boost::filesystem::path>(&conf.f_info), "Write info about samples")
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
        std::cerr << "Error: Exclude chromosomes (-x) have no effect when -b is specified. Stop" << std::endl << std::endl;
        goto print_usage_and_exit;
    }
    if (conf.mode != "count" && conf.mode != "hmm") {
        std::cerr << "Unknown mode: specify one of [count, hmm]" << std::endl << std::endl;
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
        cells[i].id = i;
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
        std::cout << "Excluding " << exclude.size() << " regions" << std::endl;
        create_fixed_bins(bins, chrom_map, conf.window, exclude, hdr);
    }
    // add last element for easy calculation of number of bins
    chrom_map.push_back((int32_t)bins.size());


    // Count in bins
    for(int i = 0; i < counts.size(); ++i)
    {
        if (!count_sorted_reads(conf.f_in[i].string(), bins, chrom_map, hdr, conf.minMapQual, counts[i], cells[i])) {
            std::cerr << "Ignore sample " << conf.f_in[i].string() << std::endl;
            counts.erase(counts.begin()+i);
            --i;
        }
    }


    // Print cell information:
    if (vm.count("info")) {
        std::cout << "Writing sample information: " << conf.f_info.string() << std::endl;
        std::ofstream out(conf.f_info.string());
        if (out.is_open()) {
            out << "sample\tcell\ttotal_reads\tduplicates\tsecondary_reads" << std::endl;
            // do not sort, so cells == counts == conf.f_in
            std::vector<CellInfo> cells2 = cells; // copy
            sort(cells2.begin(), cells2.end(), [] (CellInfo const & a, CellInfo const & b) {if (a.sample_name==b.sample_name) { return a.id < b.id;} else {return a.sample_name < b.sample_name;} } );
            for (CellInfo const & cell : cells2) {
                out << cell.sample_name << "\t";
                out << conf.f_in[cell.id].stem().string() << "\t";
                out << cell.total << "\t";
                out << cell.pcr_dups << "\t";
                out << cell.secondary << std::endl;
            }
        } else {
            std::cerr << "Cannot write to " << conf.f_info.string() << std::endl;
        }
    }


    // Just print the counts and exit
    if (conf.mode == "count") {
        std::cout << "Writing " << conf.f_out.string() << std::endl;
        std::ofstream out(conf.f_out.string());
        if (out.is_open()) {
            out << "chrom\tstart\tend\tsample\tcell\tw\tc" << std::endl;
            for(unsigned i = 0; i < counts.size(); ++i) {
                for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
                    out << hdr->target_name[bins[bin].chr];
                    out << "\t" << bins[bin].start << "\t" << bins[bin].end;
                    out << "\t" << cells[i].sample_name;
                    out << "\t" << conf.f_in[i].stem().string();
                    out << "\t" << counts[i][bin].watson_count;
                    out << "\t" << counts[i][bin].crick_count;
                    out << std::endl;
                }
            }
            out.close();
        } else {
            std::cerr << "Cannot open out file: " << conf.f_out.string() << std::endl;
            return 2;
        }
        return 0;
    }
    

    // median per sample
    for(unsigned i = 0; i < counts.size(); ++i) {
        TMedianAccumulator<unsigned int> med_acc;
        for (Counter const & count_bin : counts[i])
            med_acc(count_bin.watson_count + count_bin.crick_count);
        cells[i].median_bin_count = boost::accumulators::median(med_acc);
    }





    // find well-behaving bins
    // todo

    // estimate p
    // todo
    double p = 0.1;

    // Set up and run HMM:
    hmm::HMM<unsigned, hmm::MultiVariate<hmm::NegativeBinomial> > hmm({
        "WW", "WC", "CC", "0", "W", "C", "WWW", "WWC", "WCC", "CCC", "WWWW", "WWWC", "WWCC", "WCCC", "CCCC"});


    {
        double ini  = .03;  // Initial probability for non-diploid states
        double chng = .001; // trans prob to jump into other, non-diploid state
        double dipl = .030; // trans prob to jump back to any diploid state

        // calculated:
        double in2 = (1 - 12*ini)/3;
        double st1 = 1 - 3*dipl -11*chng;
        double st2 = 1 - 2*dipl -12*chng;
        assert(in2 > 0);
        assert(st1>0.8 && st1 <1);
        assert(st2>0.8 && st2 <1);

        //   WW    WC    CC   0     W     C     WWW   WWC   WCC   CCC  WWWW  WWWC  WWCC  WCCC  CCCC
        hmm.set_initials({
            in2,  in2,  in2,  ini,  ini,  ini,  ini,  ini,  ini,  ini,  ini,  ini,  ini,  ini,  ini});
        hmm.set_transitions({
            st2,  dipl, dipl, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, // WW
            dipl, st2,  dipl, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, // WC
            dipl, dipl, st2,  chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, // CC
            dipl, dipl, dipl, st1,  chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, // 0
            dipl, dipl, dipl, chng,  st1, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, // W
            dipl, dipl, dipl, chng, chng,  st1, chng, chng, chng, chng, chng, chng, chng, chng, chng, // C
            dipl, dipl, dipl, chng, chng, chng,  st1, chng, chng, chng, chng, chng, chng, chng, chng, // WWW
            dipl, dipl, dipl, chng, chng, chng, chng,  st1, chng, chng, chng, chng, chng, chng, chng, // WWC
            dipl, dipl, dipl, chng, chng, chng, chng, chng,  st1, chng, chng, chng, chng, chng, chng, // WCC
            dipl, dipl, dipl, chng, chng, chng, chng, chng, chng,  st1, chng, chng, chng, chng, chng, // CCC
            dipl, dipl, dipl, chng, chng, chng, chng, chng, chng, chng,  st1, chng, chng, chng, chng, // WWWW
            dipl, dipl, dipl, chng, chng, chng, chng, chng, chng, chng, chng,  st1, chng, chng, chng, // WWWC
            dipl, dipl, dipl, chng, chng, chng, chng, chng, chng, chng, chng, chng,  st1, chng, chng, // WWCC
            dipl, dipl, dipl, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng,  st1, chng, // WCCC
            dipl, dipl, dipl, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng, chng,  st1, // CCCC
        });
    }

    for (unsigned i=0; i<counts.size(); ++i)
    {
        // set n and p parameters according to a mean of sample
        double n = (double)cells[i].median_bin_count / 2 * p / (1-p);
        double z = 0.5; // mean in zero bins
        hmm.set_emissions( {\
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,2*n), hmm::NegativeBinomial(p,  z)}), // WW
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  n), hmm::NegativeBinomial(p,  n)}), // WC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  z), hmm::NegativeBinomial(p,2*n)}), // CC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  z), hmm::NegativeBinomial(p,  z)}), // 0
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  n), hmm::NegativeBinomial(p,  z)}), // W
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  z), hmm::NegativeBinomial(p,  n)}), // C
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,3*n), hmm::NegativeBinomial(p,  z)}), // WWW
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,2*n), hmm::NegativeBinomial(p,  n)}), // WWC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  n), hmm::NegativeBinomial(p,2*n)}), // WCC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  z), hmm::NegativeBinomial(p,3*n)}), // CCC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,4*n), hmm::NegativeBinomial(p,  z)}), // WWWW
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,3*n), hmm::NegativeBinomial(p,  n)}), // WWWC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,2*n), hmm::NegativeBinomial(p,2*n)}), // WWCC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  n), hmm::NegativeBinomial(p,3*n)}), // WCCC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p,  z), hmm::NegativeBinomial(p,4*n)})  // CCCC
        });
        run_HMM(hmm, counts[i], chrom_map);
    }




    // print normalized counts
    if (conf.mode == "hmm")
    {
        std::cout << "Writing file " << conf.f_out.string() << std::endl;
        std::ofstream out(conf.f_out.string());
        out << "chrom\tstart\tend\tsample\tcell\tc\tw\tclass" << std::endl;
        for(unsigned i = 0; i < counts.size(); ++i) {
            for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
                Counter & cc = counts[i][bin];
                out << hdr->target_name[bins[bin].chr];
                out << "\t" << bins[bin].start << "\t" << bins[bin].end;
                out << "\t" << cells[i].sample_name;
                out << "\t" << conf.f_in[i].stem().string();
                out << "\t" << cc.crick_count;
                out << "\t" << cc.watson_count;
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
