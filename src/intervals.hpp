#ifndef calc_bins_hpp
#define calc_bins_hpp

#include <algorithm>
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

auto interval_comp = [] (Interval const & a, Interval const & b) {
    if (a.chr == b.chr)
        // special condition: If intervals start at the same position, prefer larger one!
        if (a.start == b.start)
            return (a.end - a.start) < (b.end - b.start);
        else
            return a.start < b.start;
    else
        return a.chr < b.chr;
};


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
    // read intervals
    if (!read_exclude_file(filename, hdr, intervals))
        return false;

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
                       std::vector<Interval> const & excl,
                       bam_hdr_t* hdr)
{
    auto excl_iter = excl.begin();

    for (int32_t chrom=0; chrom<hdr->n_targets; ++chrom)
    {
        // store chrom-pointer in chrom_map
        chrom_map[chrom] = (int32_t)intervals.size();

        // skip excl. chromosomes "left" of this one
        while(excl_iter != excl.end() && excl_iter->chr < chrom)
            ++excl_iter;

        unsigned pos = 0;
        while (pos < hdr->target_len[chrom]) {

            // skip excl. bins left of pos
            while(excl_iter != excl.end() && excl_iter->chr == chrom && excl_iter->end <= pos)
                ++excl_iter;

            Interval ivl;
            ivl.chr = chrom;

            // if pos is inside an excl. interval, go to its end
            if (excl_iter != excl.end() && excl_iter->chr == chrom && pos >= excl_iter->start) {
                pos = excl_iter->end;

            } // if pos is ok but next interval is closer than binwidth
            else if (excl_iter != excl.end() && excl_iter->chr == chrom && pos+binwidth >= excl_iter->start) {
                ivl.start = pos;
                ivl.end   = std::min((int32_t)(excl_iter->start), (int32_t)(hdr->target_len[chrom]));
                intervals.push_back(ivl);
                pos = excl_iter->end;

            } // normal interval
            else {
                Interval ivl;
                ivl.chr = chrom;
                ivl.start = pos;
                ivl.end   = std::min(pos+binwidth, (unsigned)hdr->target_len[chrom]);
                intervals.push_back(ivl);
                pos += binwidth;
            }
        }
    }
    return true;
}



bool read_exclude_file(std::string const & filename, bam_hdr_t* hdr, std::vector<Interval> & intervals)
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
            if (tokIter!=tokens.end()) {
                std::string chrName = *tokIter++;
                int32_t tid = bam_name2id(hdr, chrName.c_str());
                if (tid >= 0 && tid < hdr->n_targets) {
                    Interval ivl;
                    ivl.chr = tid;
                    if (tokIter == tokens.end()) {// exclude whole chrom
                        ivl.start = 0;
                        ivl.end   = hdr->target_len[tid];
                    } else {
                        ivl.start = boost::lexical_cast<int32_t>(*tokIter++);
                        if (tokIter == tokens.end()) {
                            std::cerr << "Warning: Invalid line: " << line << std::endl;
                            continue;
                        }
                        ivl.end   = boost::lexical_cast<int32_t>(*tokIter++);
                        if (ivl.end <= ivl.start) {
                            std::cerr << "Warning: Invalid line: " << line << std::endl;
                            continue;
                        }
                    }
                    intervals.push_back(ivl);
                } else {
                    std::cerr << "Warning: Chromosome not found: " << chrName << std::endl;
                }
            }
        }
        interval_file.close();
    } else {
        std::cerr << "Error: Exclude file cannot be read: " << filename << std::endl;
        return false;
    }
    return true;
}



#endif /* calc_bins_hpp */
