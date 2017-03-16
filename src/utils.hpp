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



struct Counter {
    static const std::vector<std::string> label_names;
    static const std::map<std::string, uint8_t> label_id;
    unsigned int watson_count, crick_count;
    double watson_norm, crick_norm;
    uint8_t label;

    Counter() : watson_count(0), crick_count(0), watson_norm(0), crick_norm(0), label(0)
    {}

    bool set_label(std::string const & s) {
        auto iter = label_id.find(s);
        assert(iter != label_id.end());
        if (iter == label_id.end()) return false;
        label = iter->second;
        return true;
    }
    
    std::string get_label() const {
        assert(label < label_names.size());
        return label_names[label];
    }
};
const std::vector<std::string> Counter::label_names = {"unset", "none", "WW", "WC", "CC"};
const std::map<std::string, uint8_t> Counter::label_id = {
    {"unset",0},
    {"none", 1},
    {"WW",   2},
    {"WC",   3},
    {"CC",   4},
};




template <typename TReturn>
using TMedianAccumulator = boost::accumulators::accumulator_set<TReturn, boost::accumulators::stats<boost::accumulators::tag::median> >;


double sum(std::vector<double> const & vec)
{
    double sum = 0;
    for (double d : vec)
        sum += d;
    return(sum);
}


#endif /* utils_hpp */
