#pragma once

#include <string>
#include <unordered_map>
#include <vector> 
#include <numeric>
#include <cmath>
#include <functional>
#include <iostream> 

#include "contigs.hpp"
#include "overlaps.hpp"

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
    return (ovl.terminus_i == START && ovl.terminus_j == END) || (ovl.terminus_i == START && ovl.terminus_j == RCSTART);
}

bool is_end_match(const Overlap& ovl) {
    // Function returns true if overlap `ovl` is associated with end.
    return (ovl.terminus_i == END && ovl.terminus_j == START) || (ovl.terminus_i == END && ovl.terminus_j == RCEND);
}

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
    float total_dead_ends = 0;

    // Iterate over contigs
    int num_dead_ends;
    // Подсчет перекрытий, ассоциированных с началом
    int start_is_not_dead;
    // Подсчет перекрытий, ассоциированных с концом
    int end_is_not_dead;

    for (ContigIndex i = 0; i < contig_collection.size(); i++) {
            
        for (const auto& overlap : overlap_collection[i]) {

            start_is_not_dead += is_start_match(overlap);
            
            end_is_not_dead += is_end_match(overlap);
        }
        
        std::cout <<"rb "<< start_is_not_dead<<std::endl;
        std::cout <<"r "<< end_is_not_dead<<std::endl;

        if(start_is_not_dead!=0){start_is_not_dead=1;}
        if(end_is_not_dead!=0){end_is_not_dead=1;}

        // Calculate number of dead ends of the current contig
        num_dead_ends = num_contig_termini - start_is_not_dead - end_is_not_dead;
        // Add to `total_dead_ends`
        total_dead_ends += num_dead_ends;
        std::cout <<"rbffffff "<< total_dead_ends<<std::endl;
    }
    //std::cout <<"rbffffff "<< total_dead_ends<<std::endl;
    // Total number of termini taking account of multiplicity
    int total_termini = num_contig_termini * contig_collection.size();
    // Calculate the LQ coefficient
    float lq_coef = (1 - total_dead_ends / total_termini) * 100.0;

    return lq_coef;  // Округляем до двух знаков после запятой
}

int calc_exp_genome_size(const ContigCollection& contig_collection,
                         const OverlapCollection& overlap_collection) {
    // In this variable, total length of overlapping regions will be stored
    int total_overlap_len = 0;

    ContigIndex i = 0; 
    std::vector<Contig> contig; 

    // Get start-associated overlaps
    std::vector<Overlap> start_ovls;
    // Get end-associated overlaps
    std::vector<Overlap> end_ovls;

    // Iterate over contigs
    for (ContigIndex i = 0; i < contig_collection.size(); i++) {

        for (const Overlap& ovl : overlap_collection[i]) {
            if (is_start_match(ovl)) {
                start_ovls.push_back(ovl);
            }
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
            if (ovls.size() <= contig_collection[i].multplty) {
                // No extra overlaps. We will just add lengths of overlaps to `total_overlap_len`.
                std::copy_if(ovls.begin(), ovls.end(), std::back_inserter(ovls_to_add), not_already_counted);
            } else {
                // Some extra overlaps discovered.
                // We will consider only M longest overlaps, where M is contig's multiplicity.
                std::vector<Overlap> sorted_ovls = ovls;
                std::sort(sorted_ovls.begin(), sorted_ovls.end(), [](const Overlap& a, const Overlap& b) {
                    return a.ovl_len > b.ovl_len;
                });
                ovls_to_add.assign(sorted_ovls.begin(), sorted_ovls.begin() + contig_collection[i].multplty);
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

std::vector<Overlap> _get_start_matches(const std::vector<Overlap>& overlaps) {
    std::vector<Overlap> start_matches;
    for (const auto& overlap : overlaps) {
        if (is_start_match(overlap)) {
            start_matches.push_back(overlap);
        }
    }
    return start_matches;
}

std::vector<Overlap> _get_end_matches(const std::vector<Overlap>& overlaps) {
    std::vector<Overlap> end_matches;
    for (const auto& overlap : overlaps) {
        if (is_end_match(overlap)) {
            end_matches.push_back(overlap);
        }
    }
    return end_matches;
}

std::function<std::vector<Overlap>(const std::vector<Overlap>&)> _select_get_matches(const std::string& term) {
    if (term == "s") {
        return _get_start_matches;
    } else if (term == "e") {
        return _get_end_matches;
    } else {
        std::cerr << "Fatal error: invalid value passed to function `_get_overlaps_str_for_table` with argument `term`: `" << term << "`" << std::endl;
        std::cerr << "Please, contact the developer." << std::endl;
        exit(1);
    }
}


std::string _get_overlaps_str_for_table(const OverlapCollection& overlap_collection,
                                        const ContigCollection& contig_collection,
                                        ContigIndex key, const std::string& term) {
    // Select function for obtaining `term`-associated overlaps.
    auto get_matches = _select_get_matches(term);

    // Extract overlaps associated with `term` terminus for `key` contig
    auto overlaps = get_matches(overlap_collection[key]);

    // Convert `Overlap` instances to string representation
    if (overlaps.empty()) {
        return "-"; // no proper overlaps found
    } else {
        std::vector<std::string> match_strings; // vector for formatted strings

        for (const auto& ovl : overlaps) {
            // If contig does not match itself
            if (ovl.contig_i != ovl.contig_j) {
                // Letter for the first contig of the overlap
                std::string letter1(1, KEY2LETTER_MAP.at(ovl.terminus_i)[0]);
                std::string letter2(1, KEY2LETTER_MAP.at(ovl.terminus_j)[0]);


                // Convert and append
                match_strings.push_back("[" + letter1 + "=" + letter2 + "(" + contig_collection[ovl.contig_j].name + "); ovl=" + std::to_string(ovl.ovl_len) + "]");
            } else {
                match_strings.push_back("[Circle; ovl=" + std::to_string(ovl.ovl_len) + "]");
            }
        }

        // Separate formatted strings with spaces
        return std::accumulate(std::next(match_strings.begin()), match_strings.end(),
                               match_strings.front(),
                               [](const std::string& a, const std::string& b) {
                                   return a + " " + b;
                               });
    }
}

std::string _get_overlaps_str_for_log(const OverlapCollection& overlap_collection,
                                      const ContigCollection& contig_collection,
                                      ContigIndex key) {
    // Extract overlaps for the current contig
    const std::vector<Overlap>& overlaps = overlap_collection[key];

    if (overlaps.empty()) {
        return ""; // no proper overlaps found
    } else {
        std::string result; // string for formatted strings

        for (const Overlap& ovl : overlaps) {
            // If contig does not match itself
            if (ovl.contig_i != ovl.contig_j) {
                // Word for the first contig of the overlap
                std::string word1 = KEY2WORD_MAP.at(ovl.terminus_i);
                // Word for the second contig of the overlap
                std::string word2 = KEY2WORD_MAP.at(ovl.terminus_j);

                // Convert and append
                result += contig_collection[key].name + ": " +
                          word1 + " matches " + word2 +
                          " of " + contig_collection[ovl.contig_j].name +
                          " with overlap of " + std::to_string(ovl.ovl_len) + " bp\n";
            } else {
                // Contig is circular
                if (ovl.terminus_i == END && ovl.terminus_j == START) {
                    result += contig_collection[key].name +
                              ": contig is circular with overlap of " +
                              std::to_string(ovl.ovl_len) + " bp\n";
                }
                // Start of contig matches its own reverse-complement end
                else if (ovl.terminus_i == START && ovl.terminus_j == RCEND) {
                    result += contig_collection[key].name +
                              ": start is identical to its own rc-end with overlap of " +
                              std::to_string(ovl.ovl_len) + " bp\n";
                }
            }
        }

        return result;
    }
}

void write_adjacency_table(const ContigCollection& contig_collection,
                           const OverlapCollection& overlap_collection,
                           const std::string& outdpath) {
    // Сформировать путь к выходному файлу TSV
    std::string adj_table_fpath = outdpath + "_combinator_adjacent_contigs.tsv";

    std::cout << "Writing adjacency table to `" << adj_table_fpath << "`" << std::endl;

    // Открыть выходной файл для записи
    std::ofstream outfile(adj_table_fpath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open output file" << std::endl;
        return;
    }

    // Записать заголовок таблицы
    outfile << "#\tContig name\tLength\tCoverage\tGC(%)\tMultiplicity\tAnnotation\tStart\tEnd\n";

    std::string wrk_str;

    // Пройти по контигам и записать их свойства
    for (ContigIndex i = 0; i < contig_collection.size(); ++i) {
        const Contig& contig = contig_collection[i];

        // Порядковый номер и имя
        outfile << i + 1 << "\t" << contig.name << "\t";

        // Длина
        outfile << contig.length << "\t";

        // Покрытие
        wrk_str = (contig.cov == -1) ? "-" : std::to_string(contig.cov);
        outfile << wrk_str << "\t";

        // Содержание GC в контиге
        outfile << contig.gc_content << "\t";

        // Множество копий
        outfile << contig.multplty << "\t";

        // Пустая колонка для аннотации
        outfile << "\t";

        // Информация о найденных смежностях
        // Колонка "Start"
        wrk_str = _get_overlaps_str_for_table(overlap_collection, contig_collection, i, "s");
        outfile << wrk_str << "\t";

        // Колонка "End"
        wrk_str = _get_overlaps_str_for_table(overlap_collection, contig_collection, i, "e");
        outfile << wrk_str << "\n";
    }
}

void write_full_log(const ContigCollection& contig_collection,
                    const OverlapCollection& overlap_collection,
                    const std::string& outdpath) {
    // Сформировать путь к файлу полного журнала
    std::string log_fpath = outdpath + "_combinator_full_matching_log.txt";

    std::cout << "Writing full matching log to `" << log_fpath << "`" << std::endl;

    // Открыть файл для записи
    std::ofstream outfile(log_fpath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open output file" << std::endl;
        return;
    }

    std::string wrk_str;

    // Записать информацию о найденных смежностях
    for (ContigIndex i = 0; i < contig_collection.size(); ++i) {
        // Записать совпадения, с которыми связан начало текущего контига
        wrk_str = _get_overlaps_str_for_log(overlap_collection, contig_collection, i);
        if (!wrk_str.empty()) {
            outfile << wrk_str << "\n";
        }
    }
}
