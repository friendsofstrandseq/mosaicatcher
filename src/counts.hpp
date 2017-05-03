#ifndef counts_hpp
#define counts_hpp

#include "utils.hpp"

struct Counter {
    static const std::vector<std::string> label_names;
    static const std::map<std::string, uint8_t> label_id;
    unsigned int watson_count, crick_count;
    std::string label;
    unsigned n_supplementary;

    Counter() : watson_count(0), crick_count(0), label("NA"), n_supplementary(0)
    {}

    bool set_label(std::string const & s) {
        label = s;
        return true;
    }

    std::string get_label() const {
        return label;
    }
};


typedef std::vector<Counter> TGenomeCounts;


struct CellInfo {
    unsigned median_bin_count;
    std::string sample_name;
    int32_t id;
    // read counts
    unsigned n_mapped, n_pcr_dups, n_supplementary, n_low_mapq, n_read2s, n_counted, n_unmap;
    CellInfo() : median_bin_count(0), id(-1), n_mapped(0), n_pcr_dups(0), n_supplementary(0), n_low_mapq(0), n_read2s(0), n_counted(0), n_unmap(0) {}
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
                        TGenomeCounts & counts,
                        CellInfo & cell)
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

            if (rec->core.flag & BAM_FUNMAP) {
                continue;
            }
            ++cell.n_mapped;

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

            // Don't read every RG tag because that might slow down BAM parsing.
            // auto x = bam_aux_get(rec, "RG");

            // Ignore certain alignments
            if (rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL)) {
                ++cell.n_supplementary;
                ++(counts[bin].n_supplementary);
                continue;
            } if (rec->core.flag & BAM_FDUP) {
                ++cell.n_pcr_dups;
                continue;
            } if ((rec->core.qual < min_map_qual) || (rec->core.tid<0)) {
                ++cell.n_low_mapq;
                continue;
            } if (rec->core.flag & BAM_FREAD2) {
                ++cell.n_read2s;
                continue;
            }
            ++cell.n_counted;
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




#endif /* counts_hpp */
