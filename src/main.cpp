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
#include <htslib/sam.h>

#include "utils.hpp"

/*
    Assumptions:
    - BAM files sorted & indexed
    - All mapped to same genome
    - All belong to same sample (SM) and have distinct read groups (ID)
*/

struct Conf {
    std::vector<std::string> f_in;
    std::string f_out;
    int minMapQual = 50;
    unsigned window = 2e5;
};

struct Counter {
    int watson_count, crick_count;
    Counter() : watson_count(0), crick_count(0) {};
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
    // Test: median
    std::vector<int> test {1,2,2,2,2,2,2,4,8,20};
    std::cout << median(test) << std::endl;
    // End Test: median

    // Read arguments
    Conf c;
    c.f_out = "out.txt";
    if (argc <2) {
        std::cerr << "Please specify at least one file" << std::endl;
        return 1;
    }
    c.f_in.resize(argc - 1);
    for (int i=1; i < argc; ++i)
        c.f_in[i-1] = std::string(argv[i]);
    

    // Open all bam files
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;    
    TSamFile samfile;
    TIndex idx;
    samfile.resize(c.f_in.size());
    idx.resize(c.f_in.size());
    for(unsigned file_c = 0; file_c < c.f_in.size(); ++file_c) {
        samfile[file_c] = sam_open(c.f_in[file_c].c_str(), "r");
        if (samfile[file_c ] == NULL) {
            std::cerr << "Fail to open file " << argv[file_c] << std::endl;
            return 1;
        }
        idx[file_c] = sam_index_load(samfile[file_c], c.f_in[file_c].c_str());
        if (idx[file_c] == NULL) {
            std::cerr << "Fail to open index for " << argv[file_c] << std::endl;
            return 1;
        }
    }

    // Count in bins
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);
    typedef std::vector<std::vector<Counter> > TGenomeCounts;
    std::vector<TGenomeCounts> counts;
    counts.resize(c.f_in.size());

    for(unsigned file_c = 0; file_c < c.f_in.size(); ++file_c) {
        std::cout << "Counting: " << c.f_in[file_c] << std::endl;

        // todo: optimize. only few chroms will be used.
        counts[file_c].resize(hdr->n_targets);
        for (int chr_idx = 0; chr_idx < hdr->n_targets; ++chr_idx) {
            
            if (hdr->target_len[chr_idx] < c.window) continue;
            int bins = hdr->target_len[chr_idx] / c.window + 1;
            std::vector<Counter> & counter = counts[file_c][chr_idx];
            counter.resize(bins, Counter());

            hts_itr_t* iter = sam_itr_queryi(idx[file_c], chr_idx, 0, 
                                             hdr->target_len[chr_idx]);
            bam1_t* rec = bam_init1();

            while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
                
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

    // Print counts
    std::ofstream out(c.f_out);
    out << "chrom\tstart\tend\tsample\tw\tc" << std::endl;
    for(unsigned file_c = 0; file_c < samfile.size(); ++file_c) {
        for (int chr_idx = 0; chr_idx < hdr->n_targets; ++chr_idx) {

            if (hdr->target_len[chr_idx] < c.window) continue; 
            int bins = hdr->target_len[chr_idx] / c.window + 1;

            for (int bin = 0; bin < bins; ++bin) {
                out << hdr->target_name[chr_idx];
                out << "\t" << bin*1000000 << "\t" << (bin+1)*1000000;
                out << "\t" << c.f_in[file_c];
                out << "\t" << counts[file_c][chr_idx][bin].watson_count;
                out << "\t" << counts[file_c][chr_idx][bin].crick_count;
                out << std::endl;
            }
        }
    }
    out.close();



    // Close bam files
    for(unsigned file_c = 0; file_c < c.f_in.size(); ++file_c) {
        hts_idx_destroy(idx[file_c]);
        sam_close(samfile[file_c]);
    }


    return 0;
}