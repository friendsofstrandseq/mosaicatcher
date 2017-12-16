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
    int32_t small_intv_size;
};



int main_sces(int argc, char **argv)
{

    // Command line options
    Conf_sces conf;
    boost::program_options::options_description po_generic("Generic options");
    po_generic.add_options()
    ("help,?", "show help message")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("segments.txt"), "output file for counts")
    ("ignore-interval-size,s", boost::program_options::value<int32_t>(&conf.small_intv_size)->default_value(3000000), "Ignore segments of this size or smaller")
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
        std::cout << "Usage:   " << argv[0] << " [OPTIONS] counts.txt.gz" << std::endl << std::endl;
        std::cout << visible_options << std::endl;
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


    // cc = consecutive counts
    std::vector<std::vector<std::pair<Interval, std::string>>> cc(counts.size());
    for (unsigned i = 0; i < counts.size(); ++i) {
        for (int32_t chrom = 0; chrom < chromosomes.size(); ++chrom) {

            std::vector<std::pair<Interval, std::string>> & cci = cc[i];

            // Skip empty chromosomes
            if (chrom_map[chrom+1] - chrom_map[chrom] < 1) {
                if (vm.count("verbose") && i==0) std::cout << "[Info] Skip chromosome (" << chromosomes[chrom] << ")" << std::endl;
                continue;
            }

            // Find consecutive intervals!
            int32_t bin = chrom_map[chrom];
            cci.push_back(std::make_pair(bins[bin], counts[i][bin].label));
            for (int32_t bin = chrom_map[chrom]+1; bin < chrom_map[chrom+1]; ++bin) {
                if (counts[i][bin].label == cc[i].back().second) {
                    cci.back().first.end = bins[bin].end;
                } else {

                    cci.push_back(std::make_pair(bins[bin], counts[i][bin].label));
                }
            }


            // Get the majoriy type
            std::string majt;
            {
                unsigned ww=0, wc=0, cc=0, none=0;
                for (auto const & x : cci) {
                    if (x.second == "WW")      { ww += x.first.end - x.first.start;}
                    else if (x.second == "CC") { cc += x.first.end - x.first.start;}
                    else if (x.second == "WC") { wc += x.first.end - x.first.start;}
                    else if (x.second == "None") {none += x.first.end - x.first.start;}
                    else {
                        std::cerr << "[Warning] Unknown state " << x.second << "! Allowed states are only WW,WC,CC, and None" << std::endl;
                    }
                }
                unsigned max = std::max({ww,wc,cc});
                if (ww == max) majt = "WW";
                else if (wc == max) majt = "WC";
                else majt = "CC";
            }

            if (cci.size()==1) continue;

            // Part 1
            // Remove None and small segments.
            // At first, make sure to extend to the ends of the chromosomes.
            if (cci.size()>1 && cci[0].second == "None") {
                cci[1].first.start = cci[0].first.start;
                cci.erase(cci.begin());
            }
            unsigned N = cci.size() -1;
            if (cci.size()>1 && cci[N].second == "None") {
                cci[N-1].first.end = cci[N].first.end;
                cci.pop_back();
            }

            // Then copy by dropping small and none segments
            std::vector<std::pair<Interval, std::string>> new_cci;
            std::copy_if(cci.begin(), cci.end(),
                         std::back_inserter(new_cci),
                         [&majt, &conf](std::pair<Interval,std::string> const & x) {return x.second=="None" || (x.second != majt && x.first.end-x.first.start < conf.small_intv_size); });


            std::cout << sample_cell_names[i].second << ", " << chromosomes[chrom] << std::endl;
            for (unsigned j = 0; j < new_cci.size(); ++j) {
                Interval intv = new_cci[j].first; // copy
                std::string label = new_cci[j].second; // copy
                std::cout << "\t" << intv.start/(float)1e3 << "-" << intv.end/(float)1e3 << "\t" << label << std::endl;
            }

            // * wenn label == None
            //    * wenn rechts und links gleich, oder rechts oder links nichts mehr, dann füge die Segmente zusammen
            //    * wenn rechts und links unterschiedlich


        } // chrom
    } // i



    return 0;
}
