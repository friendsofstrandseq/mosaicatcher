//
//  calc_bins.cpp
//  strseq
//
//  Created by Sascha Meiers on 23/03/2017.
//  Copyright © 2017 Sascha Meiers. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <vector>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <htslib/sam.h>
#include "calc_bins.hpp"



struct Conf {
    boost::filesystem::path f_in;
    boost::filesystem::path f_out;
    int minMapQual;
    unsigned num_reads;
};



int main(int argc, char **argv)
{
    Conf conf;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
    ("help,?", "show help message")
    ("mapq,q", boost::program_options::value<int>(&conf.minMapQual)->default_value(10), "min mapping quality")
    ("numreads,n", boost::program_options::value<unsigned int>(&conf.num_reads)->default_value(1000), "number of reads in dynamic bin")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("bins.bed"), "output file for bins")
    ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&conf.f_in), "input bam file (WGS or merged libraries)")
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
    if (vm.count("help") || (!vm.count("input-file"))) {
        std::cout << "Usage: " << argv[0] << " [OPTIONS] <wgs.bam>" << std::endl;
        std::cout << visible_options << "\n";
        return 1;
    }

    // Prepare output file
    std::ofstream out(conf.f_out.string());
    if (!out.is_open()) {
        std::cerr << "Fail to open file " << conf.f_out.string() << std::endl;
        return 1;
    }

    // Open bam file
    samFile* samfile = sam_open(conf.f_in.string().c_str(), "r");
    if (samfile == NULL) {
        std::cerr << "Fail to open file " << conf.f_in.string() << std::endl;
        return 1;
    }
    hts_idx_t* idx = sam_index_load(samfile, conf.f_in.string().c_str());
    if (idx == NULL) {
        std::cerr << "Fail to open index for " << conf.f_in.string() << std::endl;
        return 1;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    for (int chrom = 0; chrom < hdr->n_targets; ++chrom) {

        hts_itr_t* iter = sam_itr_queryi(idx, chrom, 0, hdr->target_len[chrom]);
        bam1_t* rec = bam_init1();

        unsigned counter = 0, prev_pos = 0, prev_bin = 0;
        while (sam_itr_next(samfile, iter, rec) >= 0) {

            if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP))
                continue;
            if ((rec->core.qual < conf.minMapQual) || (rec->core.tid<0))
                continue;

            if (++counter > conf.num_reads) {
                out << hdr->target_name[chrom] << "\t" << prev_bin << "\t" << rec->core.pos << std::endl;
                prev_bin = (unsigned) rec->core.pos;
                counter = 0;
            }
            assert((unsigned)rec->core.pos >= prev_pos); // reads must be sorted
            prev_pos = (unsigned) rec->core.pos;
        }
        // end of the file: print half a bin
        if (counter >0) {
            out << hdr->target_name[chrom] << "\t" << prev_bin << "\t" << hdr->target_len[chrom] << std::endl;
        }
        bam_destroy1(rec);
        hts_itr_destroy(iter);
    }

    hts_idx_destroy(idx);
    sam_close(samfile);

    out.close();
    return 0;
}

