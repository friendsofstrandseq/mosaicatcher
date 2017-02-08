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
#include <boost/filesystem.hpp>

#include <htslib/sam.h>

#include "utils.hpp"

/*
    Assumptions:
    - BAM files sorted & indexed
    - All mapped to same genome
    - currently: one cell per file.
       - for the future: All belong to same sample (SM) and have distinct read groups (ID)
       - for the even further future: different SMs allowed.
*/

int main(int argc, char **argv)
{

    // Command line options
    Conf conf;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
        ("help,?", "show help message")
        ("mapq,q", boost::program_options::value<int>(&conf.minMapQual)->default_value(50), "min mapping quality")
        ("window,w", boost::program_options::value<unsigned int>(&conf.window)->default_value(1000000), "window size")
        ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("out.txt"), "output file for counts")
        ("mode,m", boost::program_options::value<std::string>(&conf.mode)->default_value("all"), "what to compute (raw|norm)")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || !(conf.mode == "norm" || conf.mode == "raw")) {
        std::cout << "Usage: " << argv[0] << " [OPTIONS] <strand.seq1.bam> <strand.seq2.bam> ... <strand.seqN.bam>" << std::endl;
        std::cout << visible_options << "\n";
        return 1;
    } 


    // Open all bam files
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;    
    TSamFile samfile;
    TIndex idx;
    samfile.resize(conf.f_in.size());
    idx.resize(conf.f_in.size());
    for(unsigned i = 0; i < conf.f_in.size(); ++i) {
        samfile[i] = sam_open(conf.f_in[i].string().c_str(), "r");
        if (samfile[i ] == NULL) {
            std::cerr << "Fail to open file " << conf.f_in[i].string() << std::endl;
            return 1;
        }
        idx[i] = sam_index_load(samfile[i], conf.f_in[i].string().c_str());
        if (idx[i] == NULL) {
            std::cerr << "Fail to open index for " << conf.f_in[i].string() << std::endl;
            return 1;
        }
    }

    // Count in bins
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    typedef std::vector<std::vector<Counter> > TGenomeCounts;
    std::vector<TGenomeCounts> counts;
    counts.resize(conf.f_in.size());
    for(unsigned i = 0; i < counts.size(); ++i) {
        std::cout << "Reading " << conf.f_in[i].string() << std::endl;
        counts[i].resize(hdr->n_targets);
        for (int chrom = 0; chrom < hdr->n_targets; ++chrom) {
            if (hdr->target_len[chrom] < conf.window) continue;
            int bins = hdr->target_len[chrom] / conf.window + 1;
            std::vector<Counter> & counter = counts[i][chrom];
            counter.resize(bins, Counter());

            hts_itr_t* iter = sam_itr_queryi(idx[i], chrom, 0, hdr->target_len[chrom]);
            bam1_t* rec = bam_init1();
            while (sam_itr_next(samfile[i], iter, rec) >= 0) {
                if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP))
                    continue;
                if ((rec->core.qual < conf.minMapQual) || (rec->core.tid<0)) 
                    continue;

                int32_t pos = rec->core.pos + alignmentLength(rec)/2;
                if (rec->core.flag & BAM_FREAD2) 
                    if (rec->core.flag & BAM_FREVERSE) 
                        ++( counter[(int) (pos / conf.window)].crick_count );
                    else 
                        ++( counter[(int) (pos / conf.window)].watson_count );
                // Crick = + strand, watson = - strand
                else // also for unpaired reads
                    if (rec->core.flag & BAM_FREVERSE) 
                        ++( counter[(int) (pos / conf.window)].watson_count );
                    else 
                        ++( counter[(int) (pos / conf.window)].crick_count );
            }
            bam_destroy1(rec);
            hts_itr_destroy(iter);
        }
    }
    
    // Close bam files
    for(unsigned i = 0; i < conf.f_in.size(); ++i) {
        hts_idx_destroy(idx[i]);
        sam_close(samfile[i]);
    }
 
 
    // mode "raw": print raw counts
    if (conf.mode == "raw") {
        std::ofstream out(conf.f_out.string());
        out << "chrom\tstart\tend\tsample\tw\tc" << std::endl;
        for(unsigned i = 0; i < counts.size(); ++i) {
            for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom) {
                for (unsigned bin = 0; bin < counts[0][chrom].size(); ++bin) {
                    Counter & cc = counts[0][chrom][bin];
                    out << hdr->target_name[chrom];
                    out << "\t" << bin*conf.window << "\t" << (bin+1)*conf.window;
                    out << "\t" << conf.f_in[i].string();
                    out << "\t" << cc.watson_count;
                    out << "\t" << cc.crick_count;
                    out << std::endl;
                }
            }
        }
        out.close();
        return(0);
    }
    
   
    // Normalize by sample:
    for(unsigned i = 0; i < counts.size(); ++i) {
        
        TMedianAccumulator<unsigned int> med_acc;
        for (std::vector<Counter> & count_chrom : counts[i])
            for(Counter & count_bin : count_chrom)
                med_acc(count_bin.watson_count + count_bin.crick_count);
        unsigned int sample_median = boost::accumulators::median(med_acc);
        
        if (sample_median < 20) {
            std::cout << "Ignoring " << conf.f_in[i].string() << " due to too few reads. Consider increasing window size." << std::endl;
            continue;
        }
        
        for (std::vector<Counter> & count_chrom : counts[i]) {
            for(Counter & cc : count_chrom) {
                cc.watson_norm = cc.watson_count / (double)sample_median;
                cc.crick_norm  = cc.crick_count  / (double)sample_median;
            }
        }
    }



    // To do: How to handle samples that are kicked out?

    
    // Calculate median per bin:
    std::vector<std::vector<double> > bin_medians(counts[0].size());
    for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom)
        bin_medians[chrom] = std::vector<double>(counts[0][chrom].size(), 1);
    
    if (counts.size() >= 3) {
        for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom) {
            if (counts[0][chrom].size() < 1) continue;
            for (unsigned bin = 0; bin < counts[0][chrom].size(); ++bin) {
                
                TMedianAccumulator<double> med_acc;
                for(unsigned i = 0; i < counts.size(); ++i) {
                    Counter & cc = counts[i][chrom][bin];
                    med_acc(cc.watson_norm + cc.crick_norm);
                }
                bin_medians[chrom][bin] = boost::accumulators::median(med_acc);
            }
        }
    } else {
        std::cout << "Normalizaiton per bin is skipped: too few sample" << std::endl;
    }
    
    
    // Blacklist bins with too low read count:
    for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom) {
        for (unsigned bin = 0; bin < counts[0][chrom].size(); ++bin) {
            if (bin_medians[chrom][bin] < 0.1)
                std::cout << "Black list chrom_id " << chrom << " bin " << bin << std::endl;
        }
    }
    
    // Normalize per bin
    for(unsigned i = 0; i < counts.size(); ++i) {
        for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom) {
            for (unsigned bin = 0; bin < counts[0][chrom].size(); ++bin) {
                Counter & cc = counts[i][chrom][bin];
                cc.watson_norm = cc.watson_norm / bin_medians[chrom][bin];
                cc.crick_norm  = cc.crick_norm /  bin_medians[chrom][bin];
            }
        }
    }



    // mode "norm": print normalized counts
    if (conf.mode == "norm") {
        std::ofstream out(conf.f_out.string());
        out << "chrom\tstart\tend\tsample\tw\tc" << std::endl;
        for(unsigned i = 0; i < counts.size(); ++i) {
            for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom) {
                for (unsigned bin = 0; bin < counts[0][chrom].size(); ++bin) {
                    Counter & cc = counts[0][chrom][bin];
                    out << hdr->target_name[chrom];
                    out << "\t" << bin*conf.window << "\t" << (bin+1)*conf.window;
                    out << "\t" << conf.f_in[i].string();
                    out << "\t" << cc.watson_norm;
                    out << "\t" << cc.crick_norm;
                    out << std::endl;
                }
            }
        }
        out.close();
        return(0);
    } 






    return 0;
}
