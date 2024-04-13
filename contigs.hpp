#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <unordered_map> 
#include <fstream>


std::unordered_map<char, char> _COMPL_DICT = {
    {'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'},
    {'R', 'Y'}, {'Y', 'R'}, {'S', 'S'}, {'W', 'W'},
    {'K', 'M'}, {'M', 'K'}, {'B', 'V'}, {'D', 'H'},
    {'H', 'D'}, {'V', 'B'}, {'U', 'A'}, {'N', 'N'}
};

class Contig {
public:
    // Конструктор класса Contig
    Contig(const std::string& name, int length,
           float cov, float gc_content,
           const std::string& start, const std::string& rcstart,
           const std::string& end, const std::string& rcend) :
           name(name), length(length), cov(cov),
           gc_content(gc_content), start(start),
           rcstart(rcstart), end(end), rcend(rcend), multplty(0) {}

    std::string name;     // Имя
    int length;           // Длина в bp
    float cov;            // Покрытие
    float gc_content;     // Содержание GC
    std::string start;    // Префикс длиной k этого контига
    std::string rcstart;  // Обратно-комплементарная строка start
    std::string end;      // Суффикс длиной k этого контига
    std::string rcend;    // Обратно-комплементарная строка end
    int multplty;         // Количество копий этого контига в геноме (множество)
};

typedef std::vector<Contig> ContigCollection;
typedef int ContigIndex;

std::vector<std::tuple<std::string, std::string>> fasta_generator(const std::string& filepath) {
    std::vector<std::tuple<std::string, std::string>> sequences;

    std::string curr_seq_name = "";
    std::string curr_seq = "";

    std::ifstream file(filepath);
    std::string line;

    if (file.is_open()) {
        std::getline(file, curr_seq_name);

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '>' || line[0] == '@' ) {
                // Если строка пустая или начинается с '>', это новая последовательность
                if (!curr_seq.empty()) {
                    sequences.push_back(std::make_tuple(curr_seq_name, curr_seq));
                    curr_seq.clear();
                }
                curr_seq_name = line;
            } else {
                // Иначе это часть последовательности, добавляем к текущей последовательности
                curr_seq += line;
            }
        }

        // Добавляем последнюю последовательность
        if (!curr_seq.empty()) {
            sequences.push_back(std::make_tuple(curr_seq_name, curr_seq));
        }

        file.close();
    }

    return sequences;
}

std::string format_contig_name(std::string fasta_header) {
    
    if (fasta_header[0] == '_') {
        fasta_header.erase(0, 1);
    }

    if (fasta_header[fasta_header.length() - 1] == '_') {
        fasta_header.erase(fasta_header.length() - 1, 1);
    }

    std::string contig_name = fasta_header;

    return contig_name;
}

float calc_gc_сontent(std::string sequence) {
    int gc_count = 0;

    for (char base : sequence) {
        if (base == 'G' || base == 'C' || base == 'S') {
            gc_count++;
        }
    }

    float gc_content = (static_cast<float>(gc_count) / sequence.length()) * 100;

    return gc_content;
}

std::string _get_compl_base(char base) {

    return std::string(1, _COMPL_DICT[base]);
}

std::string _rc(const std::string& seq) {
    std::string result;
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        result += _get_compl_base(*it);
    }
    return result;
}


ContigCollection get_contig_collection(const std::string& filepath, int maxk) {
    
    ContigCollection contig_collection;

    std::string curr_seq_name;
    std::string curr_seq;

    // Используем функцию fasta_generator для получения последовательных пар
    std::vector<std::tuple<std::string, std::string>> contigs = fasta_generator(filepath);

    // Выводим полученные последовательности
    for (const auto& [contig_header, contig_seq] : contigs) {

        std::string contig_name = format_contig_name(contig_header);
        float cov = 0;
        float gc_content = calc_gc_сontent(contig_seq);

        
        contig_collection.push_back(Contig(
            contig_name,
            contig_seq.length(),
            cov,
            gc_content,
            contig_seq.substr(0, maxk), // start
            _rc(contig_seq.substr(0, maxk)), // rcstart
            contig_seq.substr(contig_seq.length() - maxk, maxk), // end
            _rc(contig_seq.substr(contig_seq.length() - maxk, maxk)) // rcend
        ));
    }
    return contig_collection;
}