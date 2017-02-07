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

struct Conf {
    std::vector<boost::filesystem::path> f_in;
    boost::filesystem::path f_out;
    int minMapQual;
    unsigned int window;
    std::string mode;
};


inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
      if (bam_cigar_op(cigar[i]) == BAM_CMATCH) alen+=bam_cigar_oplen(cigar[i]);
    return alen;
  }


int main(int argc, char **argv)
{

    // Command line options
    Conf c;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
        ("help,?", "show help message")
        ("mapq,q", boost::program_options::value<int>(&c.minMapQual)->default_value(50), "min mapping quality")
        ("window,w", boost::program_options::value<unsigned int>(&c.window)->default_value(1000000), "window size")
        ("out,o", boost::program_options::value<boost::filesystem::path>(&c.f_out)->default_value("out.txt"), "output file for counts")
        ("mode,m", boost::program_options::value<std::string>(&c.mode)->default_value("all"), "what to compute (raw|norm)")
        ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file", boost::program_options::value<std::vector<boost::filesystem::path> >(&c.f_in), "input bam file(s)")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || !(c.mode == "norm" || c.mode == "raw")) {
        std::cout << "Usage: " << argv[0] << " [OPTIONS] <strand.seq1.bam> <strand.seq2.bam> ... <strand.seqN.bam>" << std::endl;
        std::cout << visible_options << "\n";
        return 1;
    } 


    // Open all bam files
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;    
    TSamFile samfile;
    TIndex idx;
    samfile.resize(c.f_in.size());
    idx.resize(c.f_in.size());
    for(unsigned i = 0; i < c.f_in.size(); ++i) {
        samfile[i] = sam_open(c.f_in[i].string().c_str(), "r");
        if (samfile[i ] == NULL) {
            std::cerr << "Fail to open file " << c.f_in[i].string() << std::endl;
            return 1;
        }
        idx[i] = sam_index_load(samfile[i], c.f_in[i].string().c_str());
        if (idx[i] == NULL) {
            std::cerr << "Fail to open index for " << c.f_in[i].string() << std::endl;
            return 1;
        }
    }

    // Count in bins
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    typedef std::vector<std::vector<Counter> > TGenomeCounts;
    std::vector<TGenomeCounts> counts;
    counts.resize(c.f_in.size());

    for(unsigned i = 0; i < c.f_in.size(); ++i) {
        std::cout << "Reading " << c.f_in[i].string() << std::endl;
        counts[i].resize(hdr->n_targets);
        for (int chr_idx = 0; chr_idx < hdr->n_targets; ++chr_idx) {


            if (hdr->target_len[chr_idx] < c.window) continue;
            int bins = hdr->target_len[chr_idx] / c.window + 1;
            std::vector<Counter> & counter = counts[i][chr_idx];
            counter.resize(bins, Counter());

            hts_itr_t* iter = sam_itr_queryi(idx[i], chr_idx, 0,
                                             hdr->target_len[chr_idx]);
            bam1_t* rec = bam_init1();

            while (sam_itr_next(samfile[i], iter, rec) >= 0) {
                if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP |\
                        BAM_FSUPPLEMENTARY | BAM_FUNMAP))
                    continue;
                if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) 
                    continue;

                int32_t pos = rec->core.pos + alignmentLength(rec)/2;

                if (rec->core.flag & BAM_FREAD1) 
                    if (rec->core.flag & BAM_FREVERSE) 
                        ++( counter[(int) (pos / c.window)].crick_count );
                    else 
                        ++( counter[(int) (pos / c.window)].watson_count );
                else
                    if (rec->core.flag & BAM_FREVERSE) 
                        ++( counter[(int) (pos / c.window)].watson_count );
                    else 
                        ++( counter[(int) (pos / c.window)].crick_count );
            }

            bam_destroy1(rec);
            hts_itr_destroy(iter);
        }
    }
    
    // Close bam files
    for(unsigned i = 0; i < c.f_in.size(); ++i) {
        hts_idx_destroy(idx[i]);
        sam_close(samfile[i]);
    }
 
 
    // mode "raw": print raw counts
    if (c.mode == "raw") {
        
        // Print counts
        std::ofstream out(c.f_out.string());
        out << "chrom\tstart\tend\tsample\tw\tc" << std::endl;
        for(unsigned i = 0; i < samfile.size(); ++i) {
            for (int chr_idx = 0; chr_idx < hdr->n_targets; ++chr_idx) {

                if (hdr->target_len[chr_idx] < c.window) continue; 
                int bins = hdr->target_len[chr_idx] / c.window + 1;

                for (int bin = 0; bin < bins; ++bin) {
                    out << hdr->target_name[chr_idx];
                    out << "\t" << bin*c.window << "\t" << (bin+1)*c.window;
                    out << "\t" << c.f_in[i].string();
                    out << "\t" << counts[i][chr_idx][bin].watson_count;
                    out << "\t" << counts[i][chr_idx][bin].crick_count;
                    out << std::endl;
                }
            }
        }
        out.close();
        return(0);
    }
    
   
    // Normalize by sample:
    for(unsigned i = 0; i < c.f_in.size(); ++i) {
        
        TMedianAccumulator<unsigned int> med_acc;
        for (std::vector<Counter> & count_chrom : counts[i])
            for(Counter & count_bin : count_chrom)
                med_acc(count_bin.watson_count + count_bin.crick_count);
        unsigned int sample_median = boost::accumulators::median(med_acc);
        
        if (sample_median < 20) {
            std::cout << "Ignoring " << c.f_in[i].string() << " due to too few reads. Consider increasing window size." << std::endl;
            continue;
        }
        
        for (std::vector<Counter> & count_chrom : counts[i]) {
            for(Counter & count_bin : count_chrom) {
                count_bin.watson_norm = count_bin.watson_count / (double)sample_median;
                count_bin.crick_norm  = count_bin.crick_count  / (double)sample_median;
            }
        }
    }

    // To do: How to handle samples that are kicked out?

    
    // Calculate median per bin:
    std::vector<std::vector<double> > bin_medians(counts[0].size());
    for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom)
        bin_medians[chrom] = std::vector<double>(counts[0][chrom].size(), 1);
    
    if (c.f_in.size() >= 3) {
        for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom) {
            if (counts[0][chrom].size() < 1) continue;
            for (unsigned bin = 0; bin < counts[0][chrom].size(); ++bin) {
                
                TMedianAccumulator<double> med_acc;
                for(unsigned i = 0; i < c.f_in.size(); ++i) {
                    Counter & c = counts[i][chrom][bin];
                    med_acc(c.watson_norm + c.crick_norm);
                }
                bin_medians[chrom][bin] = boost::accumulators::median(med_acc);
            }
        }
    } else {
        std::cout << "Normalizaiton per bin is skipped: too few sample" << std::endl;
    }
    
    
    // Blacklist bins with too low read count:
    for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom) {
        if (counts[0][chrom].size() < 1) continue;
        for (unsigned bin = 0; bin < counts[0][chrom].size(); ++bin) {
            if (bin_medians[chrom][bin] < 0.1)
                std::cout << "Black list chrom_id " << chrom << " bin " << bin << std::endl;
        }
    }
    
    // Normalize per bin
    for(unsigned i = 0; i < c.f_in.size(); ++i) {
        for (unsigned chrom = 0; chrom < counts[0].size(); ++chrom) {
            if (counts[0][chrom].size() < 1) continue;
            for (unsigned bin = 0; bin < counts[0][chrom].size(); ++bin) {
                Counter & c = counts[i][chrom][bin];
                c.watson_norm = c.watson_norm / bin_medians[chrom][bin];
                c.crick_norm  = c.crick_norm /  bin_medians[chrom][bin];
            }
        }
    }



    // mode "norm": print normalized counts
    if (c.mode == "norm") {
        
        // Print counts
        std::ofstream out(c.f_out.string());
        out << "chrom\tstart\tend\tsample\tw\tc" << std::endl;
        for(unsigned i = 0; i < samfile.size(); ++i) {
            for (int chr_idx = 0; chr_idx < hdr->n_targets; ++chr_idx) {

                if (hdr->target_len[chr_idx] < c.window) continue; 
                int bins = hdr->target_len[chr_idx] / c.window + 1;

                for (int bin = 0; bin < bins; ++bin) {
                    out << hdr->target_name[chr_idx];
                    out << "\t" << bin*c.window << "\t" << (bin+1)*c.window;
                    out << "\t" << c.f_in[i].string();
                    out << "\t" << counts[i][chr_idx][bin].watson_norm;
                    out << "\t" << counts[i][chr_idx][bin].crick_norm;
                    out << std::endl;
                }
            }
        }
        out.close();
        return(0);
    } 






    return 0;
}
