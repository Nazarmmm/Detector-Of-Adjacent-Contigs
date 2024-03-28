#pragma once

#include <vector> 
#include <algorithm>
#include <cmath>

#include "contigs.hpp"
#include "overlaps.hpp"


class CoverageCalculator {
public:
    CoverageCalculator(const ContigCollection& contig_collection) {
        _coverages = _filter_non_none_covs(contig_collection);
    }

    float get_min_coverage() const {
        float min_cov = std::numeric_limits<float>::infinity();
        for (const auto& cov : _coverages) {
            min_cov = std::min(min_cov, cov);
        }
        return std::isinf(min_cov) ? std::numeric_limits<float>::quiet_NaN() : min_cov;
    }

    float get_max_coverage() const {
        float max_cov = -std::numeric_limits<float>::infinity();
        for (const auto& cov : _coverages) {
            max_cov = std::max(max_cov, cov);
        }
        return std::isinf(max_cov) ? std::numeric_limits<float>::quiet_NaN() : max_cov;
    }

    float calc_mean_coverage() const {
        float sum = 0.0f;
        for (const auto& cov : _coverages) {
            sum += cov;
        }
        return _coverages.empty() ? std::numeric_limits<float>::quiet_NaN() : sum / _coverages.size();
    }

    float calc_median_coverage() const {
        std::vector<float> sorted_coverages = _coverages;
        std::sort(sorted_coverages.begin(), sorted_coverages.end());
        size_t size = sorted_coverages.size();
        return size % 2 == 0 ? (sorted_coverages[size / 2 - 1] + sorted_coverages[size / 2]) / 2
                             : sorted_coverages[size / 2];
    }

private:
    std::vector<float> _coverages;

    std::vector<float> _filter_non_none_covs(const ContigCollection& contig_collection) const {
        std::vector<float> coverages;
        for (const auto& contig : contig_collection) {
            if (contig.cov != std::numeric_limits<float>::quiet_NaN()) {
                coverages.push_back(contig.cov);
            }
        }
        return coverages;
    }
};

bool is_start_match(const Overlap& ovl) {
    // Function returns true if overlap `ovl` is associated with start.
    return (ovl.terminus_i == START && ovl.terminus_j == END)
           || (ovl.terminus_i == START && ovl.terminus_j == RCSTART);
}

bool is_end_match(const Overlap& ovl) {
    // Function returns true if overlap `ovl` is associated with end.
    return (ovl.terminus_i == END && ovl.terminus_j == START)
           || (ovl.terminus_i == END && ovl.terminus_j == RCEND);
}

float _calc_multiplty_by_coverage(float curr_coverage, float first_contig_coverage) {
    return curr_coverage / first_contig_coverage;
}

float _calc_multiply_by_overlaps(const std::vector<Overlap>& ovl_list) {
    // Function for calculating multiplicity of a given contig
    // based on the number of overlaps of this contig.
    // :param ovl_list: list of overlaps of current contig;

    // Count overlaps associated with start
    int num_start_matches = std::count_if(ovl_list.begin(), ovl_list.end(), is_start_match);
    // Count overlaps associated with end
    int num_end_matches = std::count_if(ovl_list.begin(), ovl_list.end(), is_end_match);

    // Obtain multiplicity based on the number of overlaps
    float multiplicity = std::max(1, std::min(num_start_matches, num_end_matches));

    return multiplicity;
}

void assign_multiplicity(ContigCollection& contig_collection, const OverlapCollection& overlap_collection) {
    // Function assigns multiplicity (copies of this contig in the genome) to contigs.
    // :param contig_collection: instance of `ContigCollection`
    //   returned by function `get_contig_collection`;
    // Function modifies `contig_collection` argument passed to it.

    // Coverage of 1st contig can be zero.
    // In this case we cannot calculate multiplicity of contigs based on coverage.
    bool first_cov_is_valid = true;

    // Set `first_cov_is_valid` to false if the first contig has no coverage.
    if (contig_collection[0].cov == 0) {                                                                                            ////////////?????????????????
        first_cov_is_valid = false;
    }

    // Set `first_cov_is_valid` to false if the first contig has zero coverage.
    // And report it.
    if (first_cov_is_valid && contig_collection[0].cov < 1e-6) {
        // Coverage of 1st contig is zero
        std::cout << "\n" << contig_collection[0].name << " has zero coverage (less than 1e-6 actually)." << std::endl;
        std::cout << "Multiplicity of contigs will be calculated based on overlaps instead of coverage.\n" << std::endl;
        first_cov_is_valid = false;
    }

    // Calculate multiplicity of contigs:
    for (size_t i = 0; i < contig_collection.size(); ++i) {
        // Validate coverage of the first contig
        bool calc_multiplicity_by_div = first_cov_is_valid;

        // Validate coverage of the current contig
        calc_multiplicity_by_div = calc_multiplicity_by_div && contig_collection[i].cov != 0;

        if (calc_multiplicity_by_div) {
            // Calculate multiplicity based on coverage
            contig_collection[i].multplty = _calc_multiplty_by_coverage(
                contig_collection[i].cov,
                contig_collection[0].cov                                                                     ////////?????????????????????
            );
        } else {
            // Calculate multiplicity based on overlaps
            contig_collection[i].multplty = _calc_multiply_by_overlaps(
                overlap_collection[i]
            );
        }
    }
}