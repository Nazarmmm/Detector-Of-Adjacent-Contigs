#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map> 

#include "contigs.hpp"


typedef int Terminus;
const Terminus START = 0;
const Terminus RCSTART = 1;
const Terminus END = 2;
const Terminus RCEND = 3;

class Overlap {
public:
    // Конструктор класса Overlap
    Overlap(ContigIndex contig_i, Terminus terminus_i,
            ContigIndex contig_j, Terminus terminus_j,
            int ovl_len) :
            contig_i(contig_i), terminus_i(terminus_i),
            contig_j(contig_j), terminus_j(terminus_j),
            ovl_len(ovl_len) {}

    // Поля класса
    ContigIndex contig_i; // индекс (ключ) первого контига
    Terminus terminus_i; // термин (первого контига) участвующий в перекрытии
    ContigIndex contig_j; // индекс (ключ) второго контига
    Terminus terminus_j; // термин (второго контига) участвующий в перекрытии
    int ovl_len; // длина перекрытия

    // Переопределение оператора преобразования в строку (аналог __repr__ в Python)
    std::string to_string() const {                                                 /////////Не используется (для тестов)
        return "<" + std::to_string(contig_i) + "-" + std::to_string(terminus_i) +
               "; " + std::to_string(contig_j) + "-" + std::to_string(terminus_j) +
               "; len=" + std::to_string(ovl_len) + ">";
    }

    // Переопределение оператора сравнения для проверки на равенство (аналог __eq__ в Python)
    bool operator==(const Overlap& other) const {
        return contig_i == other.contig_i &&
               terminus_i == other.terminus_i &&
               contig_j == other.contig_j &&
               terminus_j == other.terminus_j &&
               ovl_len == other.ovl_len;
    }

    // Переопределение оператора хеширования (аналог __hash__ в Python)
    size_t hash() const {
        size_t hash_value = 0;
        hash_combine(hash_value, contig_i);
        hash_combine(hash_value, terminus_i);
        hash_combine(hash_value, contig_j);
        hash_combine(hash_value, terminus_j);
        hash_combine(hash_value, ovl_len);
        return hash_value;
    }

private:
    // Вспомогательная функция для комбинирования хеша
    void hash_combine(size_t& seed, size_t value) const {
        seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
};

class OverlapCollection {
public:
    // Конструктор класса OverlapCollection
    OverlapCollection() {}

    // Метод для получения списка перекрытий, связанных с контигом по его ключу
    std::vector<Overlap> operator[](ContigIndex key) const {
        auto it = _collection.find(key);
        if (it != _collection.end()) {
            return it->second;
        } else {
            return std::vector<Overlap>(); // Возвращаем пустой вектор, если ключ не найден
        }
    }

    // Метод для получения размера коллекции
    size_t size() const {
        return _collection.size();
    }

    // Метод для добавления перекрытия в коллекцию
    void add_overlap(ContigIndex key, const Overlap& overlap) {
        _collection[key].push_back(overlap);
    }

    // Переопределение оператора преобразования в строку
    std::string to_string() const {                      /////////Не используется (для тестов)
        std::string result = "{";
        for (const auto& pair : _collection) {
            result += std::to_string(pair.first) + ": [";
            for (const auto& overlap : pair.second) {
                result += overlap.to_string() + ", ";
            }
            result += "], ";
        }
        result += "}";
        return result;
    }

        // Добавим методы begin и end для использования в цикле for
    auto begin() const { return std::begin(_collection); }
    auto end() const { return std::end(_collection); }

private:
    // Словарь, хранящий списки перекрытий для каждого контига
    std::unordered_map<ContigIndex, std::vector<Overlap>> _collection;
};

int find_overlap_s2s(const std::string& seq1, const std::string& seq2, int mink, int maxk) {
    // Function searches for identity between starts of seq1 and seq2.
    // Function regards overlap of length [mink, maxk].
    //
    // :param seq1: one sequence;
    // :param seq2: another sequence;
    // :param mink: minimum overlap;
    // :param maxk: maximum overlap;
    //
    // Returns 0 if overlap is less than 'mink' and
    //   length of overlap (which is <= maxk) otherwise.

    int i = mink;
    int overlap = 0;

    while (seq1.substr(0, i) == seq2.substr(0, i) && i <= maxk) {
        overlap++;
        i++;
    }
    
    return overlap == 0 ? 0 : (mink + overlap - 1);
}


int find_overlap_e2s(const std::string& seq1, const std::string& seq2, int mink, int maxk) {
    // Function searches for identity between end of seq1 and start of seq2.
    // Function regards overlap of length [mink, maxk].
    //
    // :param seq1: one sequence;
    // :param seq2: another sequence;
    // :param mink: minimum overlap;
    // :param maxk: maximum overlap;
    //
    // Returns 0 if overlap is less than 'mink' and
    //   length of overlap (which is <= maxk) otherwise.

    int ovl_len = 0;

    for (int i = mink; i <= std::min(maxk, static_cast<int>(seq1.length())); i++) {
        if (seq1.substr(seq1.length() - i) == seq2.substr(0, i)) {
            ovl_len = i;
        }
    }
    
    return ovl_len;
}


int find_overlap_e2e(const std::string& seq1, const std::string& seq2, int mink, int maxk) {
    // Function searches for identity between ends of seq1 and seq2.
    // Function regards overlap of length [mink, maxk].
    //
    // :param seq1: one sequence;
    // :param seq2: another sequence;
    // :param mink: minimum overlap;
    // :param maxk: maximum overlap;
    //
    // Returns 0 if overlap is less than 'mink' and
    //   length of overlap (which is <= maxk) otherwise.

    int i = mink;
    int overlap = 0;

    while (i <= maxk && seq1.length() >= i && seq2.length() >= i && seq1.substr(seq1.length() - i) == seq2.substr(seq2.length() - i)) {
        overlap++;
        i++;
    }

    return overlap == 0 ? 0 : (mink + overlap - 1);
}


OverlapCollection detect_adjacent_contigs(const ContigCollection& contig_collection,
                                          int mink, int maxk) {
    // Функция обнаруживает смежные контиги, сравнивая их термины.

    // Инициализация коллекции перекрытий
    OverlapCollection overlap_collection;

    // Подсчет количества контигов
    int num_contigs = contig_collection.size();

    // Итерация по контигам и сравнение их терминов с другими терминами
    
    for (ContigIndex i = 0; i < num_contigs; i++) {
        // Пропускаем контиги, длина которых меньше mink
        if (contig_collection[i].length <= mink) {
            std::cout << "\r" << i + 1 << "/" << num_contigs;
            continue;
        }

        // Сравнение начала текущего контига с концом текущего контига
        int ovl_len = find_overlap_e2s(contig_collection[i].end,
                                        contig_collection[i].start,
                                        mink, maxk);
        if (ovl_len > 0 && ovl_len < contig_collection[i].length) {
            overlap_collection.add_overlap(i, Overlap(i, END, i, START, ovl_len));
            overlap_collection.add_overlap(i, Overlap(i, START, i, END, ovl_len));
        }

        // Сравнение начала текущего контига с обратно-комплементарным концом текущего контига
        ovl_len = find_overlap_s2s(contig_collection[i].start,
                                    contig_collection[i].rcend,
                                    mink, maxk);
        if (ovl_len != 0) {
            overlap_collection.add_overlap(i, Overlap(i, START, i, RCEND, ovl_len));
            overlap_collection.add_overlap(i, Overlap(i, RCEND, i, START, ovl_len));
        }

        // Сравнение текущего контига с контигами от i+1 до N
        for (ContigIndex j = i + 1; j < num_contigs; j++) {
            // Сравнение начала i-го с концом j-го
            ovl_len = find_overlap_e2s(contig_collection[j].end,                    ///+++++++++++++
                                        contig_collection[i].start,
                                        mink, maxk);
            if (ovl_len != 0) {
                overlap_collection.add_overlap(i, Overlap(i, START, j, END, ovl_len));
                overlap_collection.add_overlap(j, Overlap(j, END, i, START, ovl_len));
            }

            // Сравнение конца i-го с началом j-го
            ovl_len = find_overlap_e2s(contig_collection[i].end,
                                        contig_collection[j].start,
                                        mink, maxk);
            if (ovl_len != 0) {
                overlap_collection.add_overlap(i, Overlap(i, END, j, START, ovl_len));
                overlap_collection.add_overlap(j, Overlap(j, START, i, END, ovl_len));
            }

            // Сравнение начала i-го с обратно-комплементарным началом j-го
            ovl_len = find_overlap_e2s(contig_collection[j].rcstart,                       //+++++++++++++++
                                        contig_collection[i].start,
                                        mink, maxk);
            if (ovl_len != 0) {
                overlap_collection.add_overlap(i, Overlap(i, START, j, RCSTART, ovl_len));
                overlap_collection.add_overlap(j, Overlap(j, START, i, RCSTART, ovl_len));
            }

            // Сравнение конца i-го с обратно-комплементарным концом j-го
            ovl_len = find_overlap_e2s(contig_collection[i].end,
                                        contig_collection[j].rcend,
                                        mink, maxk);
            if (ovl_len != 0) {
                overlap_collection.add_overlap(i, Overlap(i, END, j, RCEND, ovl_len));
                overlap_collection.add_overlap(j, Overlap(j, END, i, RCEND, ovl_len));
            }

            // Сравнение начала i-го с началом j-го
            ovl_len = find_overlap_s2s(contig_collection[i].start,                        //+++++++++++++++++++++
                                        contig_collection[j].start,
                                        mink, maxk);
            if (ovl_len != 0) {
                overlap_collection.add_overlap(i, Overlap(i, START, j, START, ovl_len));
                overlap_collection.add_overlap(j, Overlap(j, START, i, START, ovl_len));
            }

            // Сравнение конца i-го с концом j-го
            ovl_len = find_overlap_e2e(contig_collection[i].end,                          //++++++++++++++++
                                        contig_collection[j].end,
                                        mink, maxk);
            if (ovl_len != 0) {
                overlap_collection.add_overlap(i, Overlap(i, END, j, END, ovl_len));
                overlap_collection.add_overlap(j, Overlap(j, END, i, END, ovl_len));
            }

            // Сравнение начала i-го с обратно-комплементарным концом j-го
            ovl_len = find_overlap_s2s(contig_collection[i].start,                        //++++++++++++++++++++
                                        contig_collection[j].rcend,
                                        mink, maxk);
            if (ovl_len != 0) {
                overlap_collection.add_overlap(i, Overlap(i, START, j, RCEND, ovl_len));
                overlap_collection.add_overlap(j, Overlap(j, RCEND, i, START, ovl_len));
            }

            // Сравнение конца i-го с обратно-комплементарным началом j-го
            ovl_len = find_overlap_e2e(contig_collection[i].end,                         //++++++++++++++++++++++
                                        contig_collection[j].rcstart,
                                        mink, maxk);
            if (ovl_len != 0) {
                overlap_collection.add_overlap(i, Overlap(i, END, j, RCSTART, ovl_len));
                overlap_collection.add_overlap(j, Overlap(j, RCSTART, i, END, ovl_len));
            }
        }

        std::cout << "\r" << i + 1 << "/" << num_contigs;
    }

    std::cout << std::endl;

    return overlap_collection;
}