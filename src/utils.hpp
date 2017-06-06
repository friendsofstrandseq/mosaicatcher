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

#ifndef utils_hpp
#define utils_hpp


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/algorithm/string.hpp>




inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
      if (bam_cigar_op(cigar[i]) == BAM_CMATCH) alen+=bam_cigar_oplen(cigar[i]);
    return alen;
  }



struct CellInfo {
    unsigned median_bin_count;      /* raw counts */
    float mean_bin_count;           /* after removal of bad bins */
    std::string sample_name;
    int32_t id;                     /* position in conf.f_in */
    float nb_p, nb_n, nb_z;         /* NB parameters */
    bool pass_qc;                   /* set to false if no (good) SS library */
    
    /* Alignment statistics */
    unsigned n_mapped, n_pcr_dups, n_supplementary, n_low_mapq, n_read2s;
    unsigned n_counted, n_unmap;
    CellInfo() : median_bin_count(0), mean_bin_count(0), id(-1), nb_p(0),
                 nb_n(0), nb_z(0), pass_qc(true), n_mapped(0), n_pcr_dups(0),
                 n_supplementary(0), n_low_mapq(0), n_read2s(0), n_counted(0),
                 n_unmap(0)
    {}
};

struct SampleInfo {
    std::vector<float> means;
    std::vector<float> vars;
    float p;
    SampleInfo() : p(0.33) {}
};



template <typename TReturn>
using TMedianAccumulator = boost::accumulators::accumulator_set<TReturn, boost::accumulators::stats<boost::accumulators::tag::median> >;

template <typename TReturn>
using TMeanVarAccumulator = boost::accumulators::accumulator_set<TReturn, boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance> >;


double sum(std::vector<double> const & vec)
{
    double sum = 0;
    for (double d : vec)
        sum += d;
    return(sum);
}


// from Delly
inline bool get_SM_tag(std::string const& header, std::string& sample_name)
{
    std::set<std::string> smIdentifiers;
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of("\n"));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    bool rgPresent = false;
    for(;itH!=itHEnd; ++itH) {
        if (itH->find("@RG")==0) {
            TStrParts keyval;
            boost::split(keyval, *itH, boost::is_any_of("\t "));
            TStrParts::const_iterator itKV = keyval.begin();
            TStrParts::const_iterator itKVEnd = keyval.end();
            for(;itKV != itKVEnd; ++itKV) {
                size_t sp = itKV->find(":");
                if (sp != std::string::npos) {
                    std::string field = itKV->substr(0, sp);
                    if (field == "SM") {
                        rgPresent = true;
                        std::string rgSM = itKV->substr(sp+1);
                        smIdentifiers.insert(rgSM);
                    }
                }
            }
        }
    }
    if (smIdentifiers.size() == 1) {
        sample_name = *(smIdentifiers.begin());
        return true;
    } else {
        sample_name = "";
        return false;
    }
}


#endif /* utils_hpp */
