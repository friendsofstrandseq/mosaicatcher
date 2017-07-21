#include <iostream>
#include <vector>
#include <random>
#include <utility>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>

//#include <boost/math/distributions/negative_binomial.hpp>

#include "intervals.hpp"
#include "counts.hpp"
#include "iocounts.hpp"
#include "simul.hpp"


using interval::Interval;
using count::TGenomeCounts;
using count::Counter;

struct Conf {
    unsigned n_cells;
    unsigned window;
    boost::filesystem::path f_sv;
    boost::filesystem::path f_out;
    boost::filesystem::path f_sce;

    double p, min_cov, max_cov, alpha;
    unsigned sce_num;
};

template <typename TString>
bool file_exists(TString const & fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}


int main(int argc, char **argv)
{

    // Command line options
    Conf conf;
    boost::program_options::options_description po_generic("Generic options");
    po_generic.add_options()
    ("help,?", "show help message")
    ("window,w", boost::program_options::value<unsigned>(&conf.window)->default_value(200000), "window size of fixed windows")
    ("numcells,n", boost::program_options::value<unsigned>(&conf.n_cells)->default_value(10), "number of cells to simulate")
    ;

    boost::program_options::options_description po_out("Output options");
    po_out.add_options()
    ("out,o",         boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("out.txt.gz"), "output count file")
    ("sceFile,S",     boost::program_options::value<boost::filesystem::path>(&conf.f_sce), "output the positions of SCEs")
    ;


    boost::program_options::options_description po_rand("Radnomization parameters");
    po_rand.add_options()
    ("nbinom_p,p",    boost::program_options::value<double>(&conf.p)->default_value(0.8,"0.8"), "p parameter of the NB distirbution")
    ("minCoverage,c", boost::program_options::value<double>(&conf.min_cov)->default_value(10), "min. read coverage per bin")
    ("maxCoverage,C", boost::program_options::value<double>(&conf.max_cov)->default_value(60), "max. read coverage per bin")
    ("alpha,a",       boost::program_options::value<double>(&conf.alpha)->default_value(0.1,"0.1"), "noise added to all bins: mostly 0, but for a fraction alpha drawn from geometrix distribution")
    ("scesPerCell,s", boost::program_options::value<unsigned>(&conf.sce_num)->default_value(4), "Average number of SCEs per cell")
    ;


    boost::program_options::options_description po_hidden("Hidden options");
    po_hidden.add_options()
    ("sv_config_file", boost::program_options::value<boost::filesystem::path>(&conf.f_sv), "Config file for SVs (see details)")
    ;

    boost::program_options::positional_options_description po_positional;
    po_positional.add("sv_config_file", -1);

    boost::program_options::options_description po_cmdline_options;
    po_cmdline_options.add(po_generic).add(po_out).add(po_rand).add(po_hidden);
    boost::program_options::options_description po_visible_options;
    po_visible_options.add(po_generic).add(po_out).add(po_rand);
    boost::program_options::variables_map vm;

    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(po_cmdline_options).positional(po_positional).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("help") || !vm.count("sv_config_file") || !file_exists(conf.f_sv.string()))
    {
    print_usage_and_exit:
        std::cout << "Usage: " << argv[0] << " [OPTIONS] sv_config_file" << std::endl;
        std::cout << po_visible_options << std::endl;
        std::cout << std::endl;
        std::cout << "Simulate Strand-seq libraries" << std::endl;
        std::cout << "=============================" << std::endl;
        std::cout << "Simulate binned Strand-seq read counts of single cells" << std::endl;
        std::cout << "with chromosomes inherited randomly as WW,WC or CC." << std::endl;
        std::cout << "Introduce SVs from a given file and randomly add SCEs." << std::endl;
        std::cout << "To not include SVs, specify an empty file." << std::endl;
        std::cout << std::endl;
        std::cout << "The SV config file is a tab-separated file with 5 columns:" << std::endl;
        std::cout << "  - Chrom" << std::endl;
        std::cout << "  - Start" << std::endl;
        std::cout << "  - End"   << std::endl;
        std::cout << "  - SV type*" << std::endl;
        std::cout << "  - average allele frequency [0,1]" << std::endl;
        std::cout << std::endl;
        std::cout << "The allowed SV types are" << std::endl;
        std::cout << "  - het_del, hom_del" << std::endl;
        std::cout << "  - het_dup, hom_dup" << std::endl;
        std::cout << "  - het_inv, hom_inv" << std::endl;
        std::cout << "  - inv_dup" << std::endl;
        std::cout << "  - false_del (to simulate lower mappability region)" << std::endl;
        std::cout << std::endl;
        std::cout << "SV breakpoints inside a bin are modelled proportionally." << std::endl;
        std::cout << "Strand state annotation ignores SVs" << std::endl;
        return 1;
    }

    // Generate bins
    // todo: make choice of genome a parameter
    std::vector<Interval>       bins;
    std::vector<int32_t>        chrom_map;
    std::vector<int32_t>  chrom_sizes = { \
        248956422, 242193529, 198295559, 190214555,
        181538259, 170805979, 159345973, 145138636,
        138394717, 133797422, 135086622, 133275309,
        114364328, 107043718, 101991189,  90338345,
         83257441,  80373285,  58617616,  64444167,
         46709983,  50818468, 156040895,  57227415 };     // GRCh38
    std::vector<std::string> chrom_names = { \
        "chr1",  "chr2",  "chr3",   "chr4",  "chr5",  "chr6",
        "chr7",  "chr8",  "chr9",  "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
        "chr19", "chr20", "chr21", "chr22",  "chrX",  "chrY" };


    std::cout << "Resolution: " << conf.window/1000 << "kb" << std::endl;
    chrom_map.resize(chrom_sizes.size());
    create_fixed_bins(bins, chrom_map, conf.window, std::vector<Interval>(), (int32_t)chrom_sizes.size(), chrom_sizes);
    chrom_map.push_back((int32_t)bins.size());

    // Random generator
    std::random_device rd;
    std::mt19937 rd_gen(rd());

    // Global vars
    std::vector<THapCount> haplotypes;
    std::vector<THapType>  chrom_states;
    std::vector<CellInfoSimulation>  cell_info;

    // Generate basic haplotype counts of each cells, including random noise
    std::uniform_real_distribution<> rd_cov(conf.min_cov, conf.max_cov);
    std::uniform_real_distribution<> rd_unif(0,1);
    double p = conf.p;

    std::cout << "Simulating  " << conf.n_cells << " cells" << std::endl;
    for (unsigned i = 0; i < conf.n_cells; ++i)
    {
        double cov_per_bin = rd_cov(rd_gen);
        std::geometric_distribution<>         rd_geom(5/(5+log2(cov_per_bin)));
        std::negative_binomial_distribution<> rd_nb(cov_per_bin/2 * p/(1-p), p);

        CellInfoSimulation cell;
        cell.cov_per_bin = cov_per_bin;

        THapCount count(bins.size());
        for (unsigned bin = 0; bin < bins.size(); ++bin)
        {
            count[bin].h1_plus = rd_nb(rd_gen);
            count[bin].h2_plus = rd_nb(rd_gen);
            count[bin].h1_minus = (rd_unif(rd_gen)<conf.alpha) ? rd_geom(rd_gen) : 0;
            count[bin].h2_minus = (rd_unif(rd_gen)<conf.alpha) ? rd_geom(rd_gen) : 0;
        }
        haplotypes.push_back(std::move(count));
        cell_info.push_back(std::move(cell));
    }


    // SV part

    // Make sure that intervals are within bounds now!!
    std::vector<SV> sv_list;
    read_SV_config_file(conf.f_sv.string(), chrom_names, chrom_sizes, sv_list);


    // for each SV
    std::cout << "--------------------" << std::endl;
    for (SV const & sv : sv_list) {

        std::cout << "SV " << chrom_names[sv.where.chr] << ":" << sv.where.start << "-" << sv.where.end << std::endl;
        auto x = locate_partial_bins(sv.where, bins, chrom_map);

        unsigned cell_counter = 0;
        for (unsigned i = 0; i < haplotypes.size(); ++i) {

            // sample carriers
            if (rd_unif(rd_gen) < sv.vaf) {
                cell_counter++;

                float fl = x.second.first;
                float fr = x.second.second;

                if (x.first.first == x.first.second) {

                    // There is only one bin
                    flip_strand(haplotypes[i][x.first.first], sv.type, fl - (1-fr));
                } else {

                    // Partially flip first bin
                    flip_strand(haplotypes[i][x.first.first], sv.type, fl);

                    // Partially flip last bin
                    flip_strand(haplotypes[i][x.first.second], sv.type, fr);

                    // Completely flip all bins in between
                    for (unsigned bin = x.first.first+1; bin < x.first.second; ++bin) {
                        flip_strand(haplotypes[i][bin], sv.type);
                    }
                }
            }
        }
        std::cout << "   in " << cell_counter << " cells" << std::endl;
    }


    // Turn haplotypes into TGenomeCounts and simulate SCEs
    float sce_prob = (float)conf.sce_num / bins.size();
    std::vector<TGenomeCounts> final_counts;
    std::vector<unsigned> sces;
    for (unsigned i = 0; i < conf.n_cells; ++i) {
        final_counts.push_back(render_cell(haplotypes[i], chrom_map, sce_prob, sces));
    }

    // write down counts
    std::cout << "[Write] Count table " << conf.f_out.string() << std::endl;
    std::vector<std::pair<std::string, std::string>> sample_cell_names;
    for (unsigned i=0; i<conf.n_cells; ++i) {
        sample_cell_names.push_back(std::make_pair("simulated", "cell_" + std::to_string(i)));
    }
    io::write_counts_gzip(conf.f_out.string(), final_counts, bins, chrom_names, sample_cell_names);

    // write down SCEs
    if (vm.count("sceFile")) {
        std::cout << "[Write] SCE positions to file " << conf.f_sce.string() << std::endl;
        std::sort(sces.begin(), sces.end());
        for (auto x : sces)
            std::cout << bins[x] << std::endl;
    }

}
