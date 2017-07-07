#ifndef iocounts_hpp
#define iocounts_hpp

#include <iostream>
#include <vector>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "counts.hpp"       // Counter
#include "intervals.hpp"    // Interval



/** write count table in long format to a gzip file.
  *
  * @param f_out      file name to write to (should end in .gz)
  * @param counts     count matrix incl. strand state labels
  * @param bins       Lookup table for intervals, which correspond to the columns of `counts`
  * @param chrom_name Lookup table for chromosome names, e.g. `hdr->target_name`
  * @param sample_cell_name Lookup table for pair of sample and cell name, correspond to the rows of `counts`
  */
template <typename TString, typename TVec, typename TPairVec>
bool write_counts_gzip(TString const & f_out,
                       std::vector<TGenomeCounts> const & counts,
                       std::vector<Interval> const & bins,
                       TVec const & chrom_name,
                       TPairVec const & sample_cell_name)
{
    boost::iostreams::filtering_ostream out;
    boost::iostreams::file_sink fout(f_out, std::ios_base::out | std::ios_base::binary);
    out.push(boost::iostreams::gzip_compressor());
    out.push(fout);

    if (fout.is_open()) {
        out << "chrom\tstart\tend\tsample\tcell\tc\tw\tclass" << std::endl;
        for(unsigned i = 0; i < counts.size(); ++i) {
            for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
                Counter const & cc = counts[i][bin];
                out << chrom_name[bins[bin].chr];
                out << "\t" << bins[bin].start << "\t" << bins[bin].end;
                out << "\t" << sample_cell_name[i].first;
                out << "\t" << sample_cell_name[i].second;
                out << "\t" << cc.crick_count;
                out << "\t" << cc.watson_count;
                out << "\t" << cc.get_label();
                out << std::endl;
            }
        }
    } else {
        std::cerr << "[Error] Cannot open file: " << f_out << std::endl;
        return false;
    }
    return true;
}


#endif /* iocounts_hpp */
