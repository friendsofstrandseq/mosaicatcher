/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <tuple>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <htslib/sam.h>

#include "utils.hpp"
#include "version.hpp"
#include "intervals.hpp"
#include "counter.hpp"
#include "distribution.hpp"
#include "iocounts.hpp"


/**
 * @file
 * @defgroup sces SCE detection from a count table
 *
 * @todo write documentation about counting.
*/


using interval::Interval;
using count::TGenomeCounts;
using count::Counter;


struct Conf_sces {
    boost::filesystem::path f_in;
    boost::filesystem::path f_out;
    int32_t small_intv_size, low_support;
};



struct Interval2 {
    std::string label;
    unsigned start, end;
    unsigned watson_count, crick_count;
    Interval2() : label(), start(0), end(0), watson_count(0), crick_count(0) {};
    Interval2(std::string const & label, unsigned start, unsigned end, unsigned w, unsigned c) :
        label(label), start(start), end(end), watson_count(w), crick_count(c)
    {}
};


int main_sces(int argc, char **argv)
{

    // Command line options
    Conf_sces conf;
    boost::program_options::options_description po_generic("Generic options");
    po_generic.add_options()
    ("help,?", "show help message")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("segments.txt"), "output file for counts")
    ("ignore-small-regions,u", boost::program_options::value<int32_t>(&conf.small_intv_size)->implicit_value(3000000), "Ignore segments of this size or smaller")
    ("ignore-low-support-regions,v", boost::program_options::value<int32_t>(&conf.low_support)->implicit_value(100), "Ignore segments with less reads than this")
    ;

    boost::program_options::options_description po_hidden("Hidden options");
    po_hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&conf.f_in), "mosaicatcher count file (.txt.gz)")
    ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", 1);


    boost::program_options::options_description cmdline_options;
    cmdline_options.add(po_generic).add(po_hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(po_generic);
    boost::program_options::variables_map vm;

    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);


    // Check arguments
    if (vm.count("help") || !vm.count("input-file")) {

        std::cout << std::endl;
        std::cout << "Mosaicatcher " << STRINGIFYMACRO(MOSAIC_VERSION_MAJOR);
        std::cout << "." << STRINGIFYMACRO(MOSAIC_VERSION_MINOR) << std::endl;
        std::cout << "> Call Sister chromatid exchange events (SCEs)." << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:   " << argv[0] << " [-u] [-v] [OPTIONS] counts.txt.gz" << std::endl << std::endl;
        std::cout << visible_options << std::endl;
        if (!vm.count("help")) {
            std::cout << "Specify --help for more details." << std::endl;
        } else {
            std::cout <<
            "Sister Chromatid Exchange events (SCEs) are visible in Strand-seq\n"
            "data as a change in the inherited strand states (e.g. from WW to WC)\n"
            "which occurs at a random position in the chromosome. Here I implemented\n"
            "a very heuristic approach to SCE detecion (described below). It is\n"
            "recommended to provide a count table in a low resolution (e.g. 500kb)\n"
            "to reduce the number of state flips.\n"
            "\n"
            "Algorithm:\n"
            " 1. Combine consecutive intervals of the same state (taken from the\n"
            "    class column in the count table) into larger intervals.\n"
            " 2. Drop 'None' regions, which had been marked as bad regions across\n"
            "    all cells. Special care is taken at the ends of the chromosomes:\n"
            "    if a chromosome ends in a 'None' region, the adjacent region is\n"
            "    extended to the end.\n"
            " 3. Remove small (-u, recommended) or low-supported (-v) regions which\n"
            "    are not concordant with the majority type of the chromosome. This\n"
            "    should again reduce the number of state flips.\n"
            "    NOTE: Only a single type per chromosome is considered - there are\n"
            "          chromosomes with more than one SCE, which this algorithm will\n"
            "          miss."
            " 4. ...\n"
            "" << std::endl;
            std::cout << "" << std::endl;
        }
        return vm.count("help") ? 0 : 1;
    }


    /////////////////////////////////////////////////////////// global variables
    /* counts/cells */
    std::vector<std::pair<std::string,std::string>> sample_cell_names;
    std::vector<TGenomeCounts> counts;

    /* bins */
    std::vector<std::string> chromosomes;
    std::vector<Interval> bins;
    std::vector<int32_t>  chrom_map;
    /////////////////////////////////////////////////////////// global variables



    // 1. Reading count file
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::cout << "[Info] Reading file " << conf.f_in.string() << " (this is a bit slow)..." << std::endl;
    if (!io::read_counts_gzip(conf.f_in.string(),
                              counts,
                              chromosomes,
                              sample_cell_names,
                              bins))
        return 1;
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    if (vm.count("verbose")) std::cout << "[Info] " << counts.size() << " cells found." << std::endl;
    if (vm.count("verbose")) std::cout << "[Time] Reading took " << time.count() << " seconds." << std::endl;



    // Create chromosome map
    chrom_map = std::vector<int32_t>(chromosomes.size());
    if (!make_chrom_map(bins, chrom_map)) {
        std::cerr << "[Error] Something is wrong with the intervals" << std::endl;
        return 2;
    }
    chrom_map.push_back((int32_t)bins.size());
    if (vm.count("verbose")) std::cout << "[Info] " << bins.size() << " intervals across " << chromosomes.size() << " chromosomes." << std::endl;



    // 2. remove bad cells (i.e. with all None).
    std::vector<unsigned> good_cells(counts.size());
    std::iota(good_cells.begin(), good_cells.end(), 0); // 0,1,2,...
    // todo: enable this later !!
    // std::vector<unsigned> good_cells = io::get_good_cells(counts);
    unsigned prev_num_cells = counts.size();

    for (unsigned i = 0; i < good_cells.size(); ++i) {
        if (i != good_cells[i]) { // save some assignments
            counts[i] = counts[good_cells[i]];
            sample_cell_names[i] = sample_cell_names[good_cells[i]];
        }
    }
    counts.resize(good_cells.size()); // should
    std::cout << "[Warning] Removed " << prev_num_cells - counts.size() << " bad cells, leaving " << counts.size() << " good ones" << std::endl;
    if (counts.empty()) {
        std::cout << "[Error] No cells left" << std::endl;
        return 2;
    }


    // strand_states[cell][chrom][interval]
    std::vector<std::vector<std::vector<Interval2>>> strand_states(counts.size(),
                                                                   std::vector<std::vector<Interval2>>(chromosomes.size()));

    // Per sample, per chrom
    for (unsigned i = 0; i < counts.size(); ++i) {
        for (int32_t chrom = 0; chrom < chromosomes.size(); ++chrom) {


            std::vector<Interval2> cci;

            // Skip empty chromosomes
            if (chrom_map[chrom+1] - chrom_map[chrom] < 1) {
                if (vm.count("verbose") && i==0) std::cout << "[Info] Skip chromosome (" << chromosomes[chrom] << ")" << std::endl;
                continue;
            }

            // Step 1:
            // Find consecutive intervals
            int32_t bin = chrom_map[chrom];
            cci.push_back(Interval2(counts[i][bin].label,
                                   bins[bin].start,
                                   bins[bin].end,
                                   counts[i][bin].watson_count,
                                   counts[i][bin].crick_count));
            for (int32_t bin = chrom_map[chrom]+1;
                 bin < chrom_map[chrom+1]; ++bin)
            {
                if (counts[i][bin].label == cci.back().label) {
                    cci.back().end = bins[bin].end;
                    cci.back().watson_count += counts[i][bin].watson_count;
                    cci.back().crick_count += counts[i][bin].crick_count;
                } else {
                    cci.push_back(Interval2(counts[i][bin].label,
                                           bins[bin].start,
                                           bins[bin].end,
                                           counts[i][bin].watson_count,
                                           counts[i][bin].crick_count));
                }
            }

            // Step 2:
            // Get the majoriy type
            std::string majt;
            {
                unsigned ww=0, wc=0, cc=0, none=0;
                for (auto const & x : cci) {
                    if (x.label == "WW")      { ww += x.end - x.start;}
                    else if (x.label == "CC") { cc += x.end - x.start;}
                    else if (x.label == "WC") { wc += x.end - x.start;}
                    else if (x.label == "None") {none += x.end - x.start;}
                    else {
                        std::cerr << "[Warning] Unknown state " << x.label << "! Allowed states are only WW, WC, CC, and None" << std::endl;
                    }
                }
                unsigned max = std::max({ww,wc,cc});
                if (ww == max) majt = "WW";
                else if (wc == max) majt = "WC";
                else majt = "CC";
            }

            if (cci.size()==1) continue;


            // Part 3:
            // Mke sure to extend to the ends of the chromosomes by replacing
            // 'None' segments with the neighboring label. This is to make sure
            // that segments cover the whole chromosome (except if the SCE
            // occurs within an internal None segment, then there will be an
            // internal gap.
            {
                if (cci.size()>1 && cci[0].label == "None") {
                    cci[1].start = cci[0].start;
                    cci.erase(cci.begin());
                }
                unsigned N = cci.size() -1;
                if (N>0 && cci[N].label == "None") {
                    cci[N-1].end = cci[N].end;
                    cci.pop_back();
                }
            }


            // Step 4:
            // Drop none segments
            //std::vector<Interval2> new_cci;
            auto iter = std::copy_if(cci.begin(), cci.end(),
                         cci.begin(),
                         [](Interval2 const & x) { return x.label != "None"; });
            cci.erase(iter, cci.end());



            // Step 5:
            // Merge neighboring segments if they have the same state.
            iter = reduce_adjacent(cci.begin(),
                                   cci.end(),
                                   cci.begin(),
                                   [](Interval2 const & a, Interval2 const & b) -> bool {return a.label == b.label;},
                                   [](Interval2 const & a, Interval2 const & b) {Interval2 ret(a); ret.end = b.end; return ret;}
                                   );
            cci.erase(iter, cci.end());



            // Step 6:
            // Remove small non-majt regions, if they are small enough or not
            // supported by enough reads. Note that if multiple small non-majt
            // are next to another, those are merged and then evaluated.
            if (vm.count("ignore-small-regions") ||
                vm.count("ignore-low-support-regions"))
            {
                for (unsigned j = 0; j < cci.size(); ++j) {

                    Interval2 merged(cci[j]);
                    if (merged.label != majt) {

                        // if there are multiple non-majt elements; combine them
                        unsigned jj = j;
                        while(jj < cci.size() && cci[jj].label != majt) {
                            merged.end = cci[jj].end;
                            merged.label = "mixed";
                            merged.watson_count += cci[jj].watson_count;
                            merged.crick_count += cci[jj].crick_count;
                            ++jj;
                        }

                        if (vm.count("ignore-small-regions") &&
                            merged.end - merged.start <= conf.small_intv_size)
                        {
                            std::cout << " -> Remove [" << sample_cell_names[i].second << "  " << chromosomes[chrom] << ":" << merged.start/1.0e6 << "-" << merged.end/1.0e6 << " " << merged.label << "] because it's too small (" << (merged.end - merged.start)/1.0e6 << " Mb)" << std::endl;
                            for (unsigned k = j; k < jj; ++k)
                                cci[k].label = "remove";
                        }
                        if (vm.count("ignore-low-support-regions") &&
                            merged.watson_count + merged.crick_count < conf.low_support)
                        {
                            std::cout << " -> Remove [" << sample_cell_names[i].second << "  " << chromosomes[chrom] << ":" << merged.start/1.0e6 << "-" << merged.end/1.0e6 << " " << merged.label << "] because it has too few reads (" << (merged.end - merged.start)/1.0e6 << " Mb)" << std::endl;
                            for (unsigned k = j; k < jj; ++k)
                                cci[k].label = "remove";
                        }
                        j = jj;
                    }
                }

                // Now also remove the "remove"-labelled intervals...
                iter = std::remove_if(cci.begin(), cci.end(),
                                      [](Interval2 const & a) {return a.label == "remove";});
                // and then merge adjacent ones again.
                iter = reduce_adjacent(cci.begin(), iter,
                                       cci.begin(),
                                       [](Interval2 const & a, Interval2 const & b) -> bool {return a.label == b.label;},
                                       [](Interval2 const & a, Interval2 const & b) {Interval2 ret(a); ret.end = b.end; return ret;}
                                       );
                 cci.erase(iter, cci.end());
            } // vm.count("ignore-small-regions") || vm.count("ignore-low-support-regions"))


            // Save regions
            strand_states[i][chrom] = std::move(cci);

        } // chrom
    } // i


    // Step 7:
    // Recurrence !
    // todo


    std::ofstream out(conf.f_out.string());
    std::cout << "[Write] strand states (incl. SCEs): " << conf.f_out.string() << std::endl;
    if (out.is_open()) {
        out << "sample\tcell\tchrom\tstart\tend\tstate" << std::endl;
        for (unsigned i = 0; i < strand_states.size(); ++i)
            for (unsigned chrom = 0; chrom < chromosomes.size(); ++chrom)
                for (auto entry : strand_states[i][chrom])
                    out << sample_cell_names[i].first << "\t"
                        << sample_cell_names[i].second << "\t"
                        << chromosomes[chrom] << "\t"
                        << entry.start << "\t"
                        << entry.end << "\t"
                        << entry.label << std::endl;
    } else {
        std::cerr << "[Warning] Cannot write to " << conf.f_out << std::endl;
        return 2;
    }

    return 0;
}
