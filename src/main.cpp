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
#include "utils.hpp"

/*
    Assumptions:
    - BAM files sorted & indexed
    - All mapped to same genome
    - currently: one cell per file.
       - for the future: All belong to same sample (SM) and have distinct read groups (ID)
       - for the even further future: different SMs allowed.
*/

typedef std::vector<Counter> TGenomeCounts;


struct Conf {
    std::vector<boost::filesystem::path> f_in;
    boost::filesystem::path f_out;
    boost::filesystem::path f_bins;
    boost::filesystem::path f_excl;
    int minMapQual;
    unsigned int window;
    std::string mode;
};


/**
 *  count_sorted_reads
 *  ------------------
 *  Count start positions, which are expected to be sorted, into sorted bins.
 *  Suitable for both fixed and variable-width bins.
 */
bool count_sorted_reads(std::string const & filename,
                                std::vector<Interval> const & bins,
                                std::vector<int32_t> const & chrom_map,
                                bam_hdr_t * hdr,
                                int min_map_qual,
                                TGenomeCounts & counts)
{

    // Open bam file
    samFile* samfile = sam_open(filename.c_str(), "r");
    if (samfile == NULL) {
        std::cerr << "Fail to open file " << filename << std::endl;
        return false;
    }
    hts_idx_t* idx = sam_index_load(samfile, filename.c_str());
    if (idx == NULL) {
        std::cerr << "Fail to open index for " << filename << std::endl;
        return false;
    }

    std::cout << "Reading " << filename << std::endl;
    counts.resize(bins.size(), Counter());

    // access samfile chrom per chrom
    for (int32_t chrom = 0; chrom < hdr->n_targets; ++chrom) {

        // skip chromosomes with no bins
        if (chrom_map[chrom+1] - chrom_map[chrom] < 1)
            continue;

        unsigned bin = chrom_map[chrom];
        int32_t prev_pos = 0;
        hts_itr_t* iter = sam_itr_queryi(idx, chrom, 0, hdr->target_len[chrom]);
        bam1_t* rec = bam_init1();
        while (sam_itr_next(samfile, iter, rec) >= 0) {

            // Ignore certain reads
            if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP))
                continue;
            if ((rec->core.qual < min_map_qual) || (rec->core.tid<0))
                continue;
            if (rec->core.flag & BAM_FREAD2)
                continue;

            // expect pos to be sorted
            int32_t pos = rec->core.pos;
            assert(pos >= prev_pos);

            // skip all bins left of this position.
            // Stop when all bins of the chromosome are done
            while (pos >= bins[bin].end)
                if (bin++ == chrom_map[chrom+1])
                    goto end_of_chromosome;

            // Ignore reads before until we reach the start of a bin
            if (pos < bins[bin].start)
                continue;

            assert(pos >= bins[bin].start && pos < bins[bin].end);

            if (rec->core.flag & BAM_FREAD2)
                continue;
            //if (rec->core.flag & BAM_FREVERSE)
            //    ++( counter[(int) (pos / conf.window)].crick_count );
            //else
            //    ++( counter[(int) (pos / conf.window)].watson_count );
            // Crick = + strand, watson = - strand
            else // also for unpaired reads
                if (rec->core.flag & BAM_FREVERSE)
                    ++( counts[bin].watson_count );
                else
                    ++( counts[bin].crick_count );
        }
    end_of_chromosome:
        bam_destroy1(rec);
        hts_itr_destroy(iter);
    }

    hts_idx_destroy(idx);
    sam_close(samfile);
    return true;
}




int main(int argc, char **argv)
{
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
        std::cout << "  * One cell per BAM file. Sample names and read groups are ignored." << std::endl;
        std::cout << "  * For paired-end data, only read 1 is counted" << std::endl;
        return 1;
    }


    // Open first sam header
    samFile* samfile = sam_open(conf.f_in[0].string().c_str(), "r");
    if (samfile == NULL) {
        std::cerr << "Fail to open file " << conf.f_in[0].string() << std::endl;
        return 1;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);


    // Bin the genome
    std::vector<Interval> bins; // stores all intervals
    std::vector<int32_t> chrom_map(hdr->n_targets, -1);
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

    // add extra element for easier calculation of number of bins in last chromosome
    chrom_map.push_back((int32_t)bins.size());

    /*
    std::cout << "tid" << "\t" << "chr" << "\t" << "map" << "\t" << "size" << std::endl;
    for (int32_t i=0; i<hdr->n_targets; ++i)
        std::cout << i << "\t" << hdr->target_name[i] << "\t" << chrom_map[i] << "\t" << (i>0 ? chrom_map[i]-chrom_map[i-1] : -2) << std::endl;
     */


    // Count in bins: sorted pos
    std::vector<TGenomeCounts> counts;
    counts.resize(conf.f_in.size());
    for(int i = 0; i < counts.size(); ++i)
    {
        if (!count_sorted_reads(conf.f_in[i].string(), bins, chrom_map, hdr, conf.minMapQual, counts[i])) {
            std::cerr << "Ignore sample " << conf.f_in[i].string() << std::endl;
            counts.erase(counts.begin()+i);
            --i;
        }

        // todo: kick out samples if median is too low
    }

    // Print counts
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

    return(0);
}
