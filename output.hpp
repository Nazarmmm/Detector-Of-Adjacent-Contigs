#pragma once

#include <string>
#include <unordered_map>

#include "contigs.hpp"
#include "overlaps.hpp"
#include "assign_multiplicity.hpp"

// Dictionary maps `Terminus` to its "letter" representation for adjacency table.
std::unordered_map<Terminus, std::string> KEY2LETTER_MAP = {
    {START, "S"},
    {RCSTART, "rc_S"},
    {END, "E"},
    {RCEND, "rc_E"}
};

// Dictionary maps `Terminus` to its "word" representation for full log.
std::unordered_map<Terminus, std::string> KEY2WORD_MAP = {
    {START, "start"},
    {RCSTART, "rc-start"},
    {END, "end"},
    {RCEND, "rc-end"}
};

int calc_sum_contig_lengths(const ContigCollection& contig_collection) {
    int sum_length = 0;
    for (const Contig& contig : contig_collection) {
        sum_length += contig.length;
    }
    return sum_length;
}

/*
// Function to check if collection is not empty
int is_not_empty(const std::vector<Overlap>& collection) {
    return static_cast<int>(collection.size() != 0);
}*/

float calc_lq_coef(const ContigCollection& contig_collection,
                   const OverlapCollection& overlap_collection) {
    // Number of termini of a contig
    int num_contig_termini = 2;
    // Total number of dead ends taking account of multiplicity
    int total_dead_ends = 0;

    // Iterate over contigs
     for (const auto& pair : overlap_collection) {
       
       ContigIndex key = pair.first;
        const std::vector<Overlap>& overlaps = pair.second;

        // Подсчет перекрытий, ассоциированных с началом
        int start_is_not_dead = 0;
        for (const auto& overlap : overlaps) {
            if (is_start_match(overlap)) {
                start_is_not_dead++;
            }
        }

        // Подсчет перекрытий, ассоциированных с концом
        int end_is_not_dead = 0;
        for (const auto& overlap : overlaps) {
            if (is_end_match(overlap)) {
                start_is_not_dead++;
            }
        }


        // Calculate number of dead ends of the current contig
        int num_dead_ends = num_contig_termini - start_is_not_dead - end_is_not_dead;
        // Add to `total_dead_ends`
        total_dead_ends += num_dead_ends;
    }

    // Total number of termini taking account of multiplicity
    int total_termini = num_contig_termini * contig_collection.size();
    // Calculate the LQ coefficient
    float lq_coef = (1 - static_cast<float>(total_dead_ends) / total_termini) * 100;

    return lq_coef;  // Округляем до двух знаков после запятой
}

int calc_exp_genome_size(const ContigCollection& contig_collection,
                         const OverlapCollection& overlap_collection) {
    // In this variable, total length of overlapping regions will be stored
    int total_overlap_len = 0;

    ContigIndex i; 
    std::vector<Contig> contig; 

    // Iterate over contigs
    for (const auto& contig : contig_collection) {

        // Get start-associated overlaps
        std::vector<Overlap> start_ovls;
        for (const Overlap& ovl : overlap_collection[i]) {
            if (is_start_match(ovl)) {
                start_ovls.push_back(ovl);
            }
        }


        // Get end-associated overlaps
        std::vector<Overlap> end_ovls;
        for (const Overlap& ovl : overlap_collection[i]) {
            if (is_end_match(ovl)) {
                end_ovls.push_back(ovl);
            }
        }

        // Function for `filter`.
        // Purpose: we won't consider contigs with index > i
        //   in order not to count an overlap twice.
        auto not_already_counted = [&](const Overlap& x) { return x.contig_j >= i; };

        for (const auto& ovls : {start_ovls, end_ovls}) {
            std::vector<Overlap> ovls_to_add;
            if (ovls.size() <= contig.multplty) {
                // No extra overlaps. We will just add lengths of overlaps to `total_overlap_len`.
                std::copy_if(ovls.begin(), ovls.end(), std::back_inserter(ovls_to_add), not_already_counted);
            } else {
                // Some extra overlaps discovered.
                // We will consider only M longest overlaps, where M is contig's multiplicity.
                std::vector<Overlap> sorted_ovls = ovls;
                std::sort(sorted_ovls.begin(), sorted_ovls.end(), [](const Overlap& a, const Overlap& b) {
                    return a.ovl_len > b.ovl_len;
                });
                ovls_to_add.assign(sorted_ovls.begin(), sorted_ovls.begin() + contig.multplty);
                std::copy_if(ovls_to_add.begin(), ovls_to_add.end(), std::back_inserter(ovls_to_add), not_already_counted);
            }

            // Add lengths of overlaps to `total_overlap_len`
            for (const auto& ovl : ovls_to_add) {
                total_overlap_len += ovl.ovl_len;
            }
        }
    }

    // Calculate length of the genome, taking account of multiplicity of contigs.
    int expected_genome_size = 0;
    for (const auto& contig : contig_collection) {
        expected_genome_size += contig.length * contig.multplty;
    }
    expected_genome_size -= total_overlap_len; // Subtract total length of overlapping regions

    return expected_genome_size;
}

void write_summary(const ContigCollection& contig_collection, const OverlapCollection& overlap_collection,
                   const std::string& infpath, const std::string& outdpath) {
    // Путь к файлу сводки
    std::string summary_fpath = outdpath;
    std::cout << "Writing summary to `" << summary_fpath << "`" << std::endl;

    // Открытие файла на запись
    std::ofstream outfile(summary_fpath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file for writing: " << summary_fpath << std::endl;
        return;
    }

    // Запись пути к входному файлу
    outfile << "Input file: `" << infpath << "`\n\n";

    // Запись сводки с некоторыми статистическими данными
    outfile << " === Summary ===\n";

    outfile << contig_collection.size() << " contigs were processed.\n";

    outfile << "Sum of contig lengths: " << calc_sum_contig_lengths(contig_collection) << " bp\n";

    outfile << "Expected length of the genome: " << calc_exp_genome_size(contig_collection, overlap_collection) << " bp\n";
    
    CoverageCalculator cov_calc(contig_collection);

    // Min coverage
    float min_coverage = cov_calc.get_min_coverage();
    outfile << "Min coverage: " << std::to_string(min_coverage) << "\n";

    // Max coverage
    float max_coverage = cov_calc.get_max_coverage();
    outfile << "Max coverage: " << std::to_string(max_coverage) << "\n";

    // Mean coverage
    float mean_coverage = cov_calc.calc_mean_coverage();
    outfile << "Mean coverage: " << std::to_string(mean_coverage) << "\n";
    
    // Median coverage
    float median_coverage = cov_calc.calc_median_coverage();
    outfile << "Median coverage: " << std::to_string(median_coverage) << "\n";
    
    outfile << "LQ-coefficient: " << calc_lq_coef(contig_collection, overlap_collection) << "\n";

}