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

enum SV_type {
    het_inv,
    hom_inv,
    het_del,
    hom_del,
    het_dup,
    hom_dup,
    inv_dup,
    false_del
};

struct SV {
    Interval where;
    SV_type type;
    float vaf;
    SV(Interval const & intvl, SV_type const & type) :
    where(intvl), type(type), vaf(1)
    {}
    SV(Interval const & intvl, SV_type const & type, float vaf) :
    where(intvl), type(type), vaf(vaf)
    {}
};


struct CellInfoSimulation {
    double cov_per_bin;
};



typedef std::vector<Hapl> THapCount;
typedef std::vector<std::string> THapType;



/** Pick a random chromosome, weighted by size of the chromosome
  *
  * @param chrom_map sorted list of indices pointing to "bins". chrom_map[3] = 8
  *        means that the first bin of chromosome "3" is bins[8].
  * @param rd generator such as std::mt19937
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




/** (Partially) flip haplotype counts according to a certail SV type
  *
  * For a het inversion for example, plus and minus on h1 are exchanged.
  * If a breakpoint of the SV does not perfectly align with bin boundaries,
  * you can specify an additional fraction f [0,1] to say that the flip should
  * only occur in f% of the bin.
  *
  * @param h Haplotype of a single bin (containing 4 counts).
  * @param sv_type Type of SV, only a few are possible.
  * @param f Apply the flip only to a portion of this bin.
  */
inline void flip_strand(Hapl & h, SV_type sv_type, float f = 1)
{
    Hapl x = h;
    switch(sv_type) {
        case het_inv:
            h.h1_plus  = (1-f) * h.h1_plus  + f * x.h1_minus;
            h.h1_minus = (1-f) * h.h1_minus + f * x.h1_plus;
            break;
        case hom_inv:
            h.h1_plus  = (1-f) * h.h1_plus  + f * x.h1_minus;
            h.h1_minus = (1-f) * h.h1_minus + f * x.h1_plus;
            h.h2_plus  = (1-f) * h.h2_plus  + f * x.h2_minus;
            h.h2_minus = (1-f) * h.h2_minus + f * x.h2_plus;
            break;
        case het_del:
            h.h1_plus  = (1-f) * h.h1_plus  + f * 0;
            h.h1_minus = (1-f) * h.h1_minus + f * 0;
            break;
        case hom_del:
            h.h1_plus  = (1-f) * h.h1_plus  + f * 0;
            h.h1_minus = (1-f) * h.h1_minus + f * 0;
            h.h2_plus  = (1-f) * h.h2_plus  + f * 0;
            h.h2_minus = (1-f) * h.h2_minus + f * 0;
            break;
        case het_dup:
            h.h1_plus  = (1-f) * h.h1_plus  + f * 2 * h.h1_plus;
            h.h1_minus = (1-f) * h.h1_minus + f * 2 * h.h1_minus; // should be 0
            break;
        case hom_dup:
            h.h1_plus  = (1-f) * h.h1_plus  + f * 2 * h.h1_plus;
            h.h1_minus = (1-f) * h.h1_minus + f * 2 * h.h1_minus; // should be 0
            h.h2_plus  = (1-f) * h.h2_plus  + f * 2 * h.h2_plus;
            h.h2_minus = (1-f) * h.h2_minus + f * 2 * h.h2_minus; // should be 0
            break;
        case inv_dup:
            h.h1_plus  = (1-f) * h.h1_plus  + f * (x.h1_minus + x.h1_plus);
            h.h1_minus = (1-f) * h.h1_minus + f * (x.h1_minus + x.h1_plus);
            break;
        case false_del:
            h.h1_plus  = (1-f) * h.h1_plus  + f * x.h1_plus/2;
            h.h1_minus = (1-f) * h.h1_minus + f * x.h1_minus/2;
            h.h2_plus  = (1-f) * h.h2_plus  + f * x.h2_plus/2;
            h.h2_minus = (1-f) * h.h2_minus + f * x.h2_minus/2;
    }
}




/** Given genomic coordinates, find the start and end bins (inclusive)
  *
  * such that bin.start <= pos <= bin.end is true for both start and end 
  * coordinate of `where`.
  *
  * @param where Interval of SV.
  * @param bins Chromosmal bins (sorted).
  * @param chrom_map Mapping to first bin of each chromosome.
  */
inline std::pair<int32_t, int32_t> locate_bins(Interval const & where,
                                    std::vector<Interval> const & bins,
                                    std::vector<int32_t> const & chrom_map)
{
    // Get bins
    float start = std::upper_bound(bins.begin(), bins.end(), Interval(where.chr, where.start, where.start)) - bins.begin() - 1;
    float end   = std::lower_bound(bins.begin(), bins.end(), Interval(where.chr, where.end,   where.end)) - bins.begin() - 1;

    assert(bins[start].chr == where.chr);
    assert(bins[start].start <= where.start && bins[start].end >= where.start);
    assert(bins[end].chr == where.chr);
    assert(bins[end].end >= where.end && bins[end].end >= where.start);

    return std::make_pair(start, end);
}

inline float left_frac(Interval const & bin, int32_t pos) {
    assert(pos >= bin.start);
    assert(pos <= bin.end);
    assert(bin.end > bin.start);
    return (float)(pos - bin.start) / (bin.end - bin.start);
}

inline float right_frac(Interval const & bin, int32_t pos) {
    assert(pos >= bin.start);
    assert(pos <= bin.end);
    assert(bin.end > bin.start);
    return (float)(bin.end - pos) / (bin.end - bin.start);
}


inline std::pair<std::pair<int32_t,int32_t>, std::pair<float,float>> locate_partial_bins(
    Interval const & where,
    std::vector<Interval> const & bins,
    std::vector<int32_t> const & chrom_map)
{
    // Get bins
    std::pair<std::pair<int32_t,int32_t>, std::pair<float,float>> tuple;
    tuple.first = locate_bins(where, bins, chrom_map);
    // determine the portion of the bins that are within the SV:
    // right_frac of start bin, and left_frac of end bin
    tuple.second = std::make_pair(right_frac(bins[tuple.first.first], where.start),
                                  left_frac(bins[tuple.first.second], where.end));
    return tuple;
}





/** Read the SV config file (5 columns)
  *
  * this function also checks that the SV definition is valid and within the 
  * chromosome boarders. Note that this is important because this check will not
  * be done during `locate_partial_bins`.
  *
  * @param filename file name.
  * @param chrom_names vector of names of the chromosomes.
  * @param chrom_size  vector of chromosome sizes. Must match `chrom_names` in size.
  * @param sv_list List of SVs to be written (push_back is used)
  */
bool read_SV_config_file(std::string const & filename,
                         std::vector<std::string> const & chrom_names,
                         std::vector<int32_t> const & chrom_size,
                         std::vector<SV> & sv_list)
{
    std::ifstream interval_file(filename.c_str(), std::ifstream::in);
    if (interval_file.is_open()) {
        while (interval_file.good()) {
            std::string line;
            getline(interval_file, line);
            typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
            boost::char_separator<char> sep(" \t,;");
            Tokenizer tokens(line, sep);
            Tokenizer::iterator tokIter = tokens.begin();
            if (tokIter!=tokens.end())
            {
                std::string chrName = *tokIter++;

                // get chromosome id
                //int32_t tid = bam_name2id(hdr, chrName.c_str());
                int32_t tid = -1;
                for (int32_t i = 0; i < chrom_names.size(); ++i) {
                    if (chrom_names[i] == chrName) {
                        tid = i;
                        break;
                    }
                }
                if (tid >= 0 && tid < chrom_names.size())
                {
                    Interval ivl;
                    ivl.chr = tid;

                    if (tokIter == tokens.end()) {
                        std::cerr << "Warning: Invalid line: " << line << std::endl;
                        continue;
                    }
                    ivl.start = boost::lexical_cast<int32_t>(*tokIter++);
                    if (tokIter == tokens.end()) {
                        std::cerr << "Warning: Invalid line: " << line << std::endl;
                        continue;
                    }
                    ivl.end   = boost::lexical_cast<int32_t>(*tokIter++);

                    // check whether interval makes sens
                    if (ivl.start < 0 || ivl.start > ivl.end || ivl.end > chrom_size[tid]) {
                        std::cerr << "[Warning] Interval out of bounds " << line << std::endl;
                        continue;
                    }

                    // Get SV type and VAF
                    if (tokIter == tokens.end()) {
                        std::cerr << "Warning: Invalid line: " << line << std::endl;
                        continue;
                    }
                    std::string type_in = *tokIter++;
                    SV_type type_out;
                    if      (type_in == "het_inv")
                        type_out = het_inv;
                    else if (type_in == "hom_inv")
                        type_out = hom_inv;
                    else if (type_in == "het_del")
                        type_out = het_del;
                    else if (type_in == "hom_del")
                        type_out = hom_del;
                    else if (type_in == "het_dup")
                        type_out = het_dup;
                    else if (type_in == "hom_dup")
                        type_out = hom_dup;
                    else if (type_in == "inv_dup")
                        type_out = inv_dup;
                    else if (type_in == "false_del")
                        type_out = false_del;
                    else {
                        std::cerr << "[Warning] Unknown SV type. Ignored: " << line << std::endl;
                    }


                    // Variant allele frequency
                    if (tokIter == tokens.end()) {
                        std::cerr << "Warning: Invalid line: " << line << std::endl;
                        continue;
                    }
                    float sv_vaf = boost::lexical_cast<float>(*tokIter++);
                    assert(tokIter == tokens.end());

                    sv_list.push_back(SV(ivl, type_out, sv_vaf));
                } else {
                    std::cerr << "Warning: Chromosome not found: " << chrName << " in \"" << line << "\"" << std::endl;
                }
            }
        }
        interval_file.close();
    } else {
        std::cerr << "[Error] SV config file cannot be read: " << filename << std::endl;
        return false;
    }
    return true;
}





#endif /* simul_hpp */
