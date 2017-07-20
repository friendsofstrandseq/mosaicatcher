#ifndef simul_hpp
#define simul_hpp

#include <iostream>
#include <vector>
#include <random>

#include "intervals.hpp"
#include "counts.hpp"

using interval::Interval;
using count::TGenomeCounts;
using count::Counter;



struct Hapl {
    unsigned h1_plus, h1_minus, h2_plus, h2_minus;
    Hapl() : h1_plus(0), h1_minus(0), h2_plus(0), h2_minus(0)
    {}
};
typedef std::vector<Hapl> THapCount;
typedef std::vector<std::string> THapType;


struct CellInfoSimulation {
    double cov_per_bin;
};


/** Pick a random chromosome, weighted by size of the chromosome
  *
  * @param chrom_map sorted list of indices pointing to "bins". chrom_map[3] = 8
  *        means that the first bin of chromosome "3" is bins[8].
  * @param random generator such as std::mt19937
  */
template <typename TRandomGenerator>
int32_t pick_random_chrom(std::vector<int32_t> chrom_map, TRandomGenerator rd)
{
    std::uniform_int_distribution<> rd_chrom(0, chrom_map.back()-1);
    int32_t rand_bin = rd_chrom(rd);
    for (int32_t i = 0; i < chrom_map.size(); ++i)
        if (chrom_map[i] >= rand_bin)
            return i;
    return -1;
}



/** Turn haplotype information into GenomeCounts
  *
  * @param hapls table of haplotypes (THapl)
  */
TGenomeCounts render_cell(THapCount const & hapls,
                          std::vector<int32_t> const & chrom_map,
                          float sce_prob,
                          std::vector<unsigned> & sce_positions)
{
    std::random_device rd;
    std::mt19937 rd_gen(rd());
    std::uniform_real_distribution<> rd_unif(0,1);

    TGenomeCounts counts(chrom_map.back());

    // strand state: 0 = WW, 1 == WC, 2 == CC
    int state = 1;
    int32_t chrom = 0;
    for (unsigned bin = 0; bin < counts.size(); ++bin)
    {
        // at the start of a chromosome: randomly select a strand state
        if (bin == (unsigned)chrom_map[chrom]) {
            state = rd_unif(rd_gen) < 0.5 ? (rd_unif(rd_gen) < 0.5 ? 0 : 2) : 1;
            chrom++;
        }

        // If an SCE occurs, change the state to one of the possible states
        if (rd_unif(rd_gen) < sce_prob) {
            if (state == 0 || state == 2) state = 1;
            else state = rd_unif(rd_gen) < 0.5 ? 0 : 2;
            sce_positions.push_back(bin);
        }

        // Now fill bins according to that state
        if (state == 0) {
            counts[bin].label = "WW";
            counts[bin].watson_count = hapls[bin].h1_plus + hapls[bin].h2_plus;
            counts[bin].crick_count  = hapls[bin].h1_minus + hapls[bin].h2_minus;
        }
        if (state == 1) {
            counts[bin].label = "WC";
            counts[bin].watson_count = hapls[bin].h1_plus + hapls[bin].h2_minus;
            counts[bin].crick_count  = hapls[bin].h2_plus + hapls[bin].h1_minus;
        }
        if (state == 2) {
            counts[bin].label = "CC";
            counts[bin].watson_count = hapls[bin].h1_minus + hapls[bin].h2_minus;
            counts[bin].crick_count  = hapls[bin].h1_plus + hapls[bin].h2_plus;
        }
    }
    return counts;
}


inline bool simulate_het_inv(THapCount & hap, unsigned bin, unsigned len) {
    for (unsigned i = bin; i < bin+len; ++i) {
        unsigned h1_plus = hap[i].h1_plus;
        hap[i].h1_plus = 0;
        hap[i].h1_minus = h1_plus;
    }
    return true;
}

inline bool simulate_hom_inv(THapCount & hap, unsigned bin, unsigned len) {
    for (unsigned i = bin; i < bin+len; ++i) {
        unsigned h1_plus = hap[i].h1_plus;
        unsigned h2_plus = hap[i].h2_plus;
        hap[i].h1_plus = 0;
        hap[i].h2_plus = 0;
        hap[i].h1_minus = h1_plus;
        hap[i].h2_minus = h2_plus;
    }
    return true;
}

inline bool simulate_het_del(THapCount & hap, unsigned bin, unsigned len) {
    for (unsigned i = bin; i < bin+len; ++i) {
        hap[i].h1_plus = 0;
    }
    return true;
}

inline bool simulate_hom_del(THapCount & hap, unsigned bin, unsigned len) {
    for (unsigned i = bin; i < bin+len; ++i) {
        hap[i].h1_plus = 0;
        hap[i].h2_plus = 0;
    }
    return true;
}

inline bool simulate_het_dup(THapCount & hap, unsigned bin, unsigned len) {
    for (unsigned i = bin; i < bin+len; ++i) {
        hap[i].h1_plus = 2 * hap[i].h1_plus;
    }
    return true;
}

inline bool simulate_hom_dup(THapCount & hap, unsigned bin, unsigned len) {
    for (unsigned i = bin; i < bin+len; ++i) {
        hap[i].h1_plus = 2 * hap[i].h1_plus;
        hap[i].h2_plus = 2 * hap[i].h2_plus;
    }
    return true;
}



struct SV {
    Interval where;
    std::string type;
    float vaf;
    SV(Interval const & intvl, std::string const & type) :
        where(intvl), type(type), vaf(1)
    {}
    SV(Interval const & intvl, std::string const & type, float vaf) :
        where(intvl), type(type), vaf(vaf)
    {}
};

std::pair<float, float> locate_bins(SV const & sv, std::vector<Interval> const & bins,
                                         std::vector<int32_t> const & chrom_map)
{
    float start = bins.begin() - std::lower_bound(bins.begin(), bins.end(), Interval(sv.where.chr, sv.where.start, sv.where.start));
    float end   = bins.begin() - std::upper_bound(bins.begin(), bins.end(), Interval(sv.where.chr, sv.where.end,   sv.where.end));

    return std::make_pair(start, end);
}





#endif /* simul_hpp */
