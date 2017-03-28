#ifndef calc_bins_hpp
#define calc_bins_hpp

#include <boost/tokenizer.hpp>
#include <htslib/sam.h>

// BED interval: 0-based, right-exclusive
struct Interval {
    int32_t chr;
    int32_t start;
    int32_t end;
    Interval() : chr(0), start(0), end(0) {}
    Interval(int32_t tid, int32_t s, int32_t e): chr(tid), start(s), end(e) {}
};

auto interval_comp = [] (Interval const & a, Interval const & b) { return (a.chr==b.chr) ? a.start < b.start : a.chr < b.chr; };


/**
 *  read_dynamic_bins
 *  -----------------
 *  read bed file. Sort intervals by chrom, start. Check they don't overlap.
 *  return sorted bin vector plus a chrom_map that stores first bin number 
 *  of each tid
 */
template <typename TFilename>
bool read_dynamic_bins(std::vector<Interval> & intervals,
                       std::vector<int32_t> & chrom_map,
                       TFilename const & filename,
                       bam_hdr_t* hdr)
{
    std::ifstream interval_file(filename, std::ifstream::in);
    if (interval_file.is_open()) {
        while (interval_file.good()) {
            std::string intervalLine;
            getline(interval_file, intervalLine);
            typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
            boost::char_separator<char> sep(" \t,;");
            Tokenizer tokens(intervalLine, sep);
            Tokenizer::iterator tokIter = tokens.begin();
            if (tokIter!=tokens.end()) {
                std::string chrName=*tokIter++;
                int32_t tid = bam_name2id(hdr, chrName.c_str());
                if (tid >= 0) {
                    if (tokIter!=tokens.end()) {
                        Interval bed;
                        bed.chr = tid;
                        bed.start = boost::lexical_cast<int32_t>(*tokIter++);
                        bed.end = boost::lexical_cast<int32_t>(*tokIter++);
                        intervals.push_back(bed);
                    }
                } else {
                    std::cerr << "chromosome not found: " << chrName << " (ignored)" << std::endl;
                }
            }
        }
        interval_file.close();
    } else {
        std::cerr << "Cannot open file " << filename << std::endl;
        return false;
    }
    if (intervals.size()<1) {
        std::cerr << "No intervals" << std::endl;
        return false;
    }

    // sort intervals
    std::sort(intervals.begin(), intervals.end(), interval_comp);

    // check that intervals don't overlap
    int32_t prev = -1, j = 0;
    for (unsigned i=0; i<intervals.size();++i) {
        if (intervals[i].chr == prev) {
            if (intervals[i-1].end > intervals[i].start) {
                std::cerr << "Intervals overlap. This is not supported!" << std::endl;
                return false;
            }
        } else {
            while (j<=intervals[i].chr)
                chrom_map[j++] = i;
            prev = intervals[i].chr;
        }
    }
    while (j<hdr->n_targets)
        chrom_map[j++] = (int32_t)intervals.size();

    return true;
}


bool create_fixed_bins(std::vector<Interval> & intervals,
                       std::vector<int32_t> & chrom_map,
                       unsigned binwidth,
                       std::vector<int32_t> const & tids, // sorted
                       bam_hdr_t* hdr)
{
    if (tids.size()<1) {
        std::cerr << "No chromosomes left" << std::endl;
        return false;
    }
    int32_t j=0;
    for (int32_t chrom : tids) {
        assert(chrom >= 0 && chrom < hdr->n_targets);

        // store chrom-pointer in chrom_map
        while (j<=chrom)
            chrom_map[j++] = (int32_t)intervals.size();
        
        // fill intervals with bins for this chromosome
        int32_t prev_pos = 0;
        for (int32_t pos = binwidth; pos < hdr->target_len[chrom]; pos += binwidth) {
            intervals.push_back(Interval(chrom, prev_pos, pos));
            prev_pos = pos;
        }
        if (prev_pos + binwidth < hdr->target_len[chrom])
            intervals.push_back(Interval(chrom, prev_pos+binwidth, hdr->target_len[chrom]));
    }
    while (j < hdr->n_targets)
        chrom_map[j++] = (int32_t)intervals.size();
    return true;
}

#endif /* calc_bins_hpp */
