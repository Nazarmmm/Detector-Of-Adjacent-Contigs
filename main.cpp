#include <iostream> //fasta bondage
#include <string>
#include <chrono>

#include "contigs.hpp"                                //2 files cra, c,сборщик spades оставляет концы равное к-меру, fastg bondage
#include "overlaps.hpp"
#include "assign_multiplicity.hpp"
#include "output.hpp"


int main() {
    // Зафиксируем начальное время
    auto start_time = std::chrono::high_resolution_clock::now();

    std::string filepath = "testg.fasta";
    int maxk=100;
    int mink=5;
    
    // Получение коллекции контигов
    ContigCollection contig_collection = get_contig_collection(filepath, maxk);

    /*for (const auto& contig : contig_collection) {
        std::cout << "Contig Name: " << contig.name << std::endl;
        std::cout << "Length: " << contig.length << std::endl;
        std::cout << "Coverage: " << contig.cov << std::endl;
        std::cout << "GC Content: " << contig.gc_content << std::endl;
        std::cout << "Start: " << contig.start << std::endl;
        std::cout << "Reverse-Complement of Start: " << contig.rcstart << std::endl;
        std::cout << "End: " << contig.end << std::endl;
        std::cout << "Reverse-Complement of End: " << contig.rcend << std::endl;
        std::cout << "Multiplicity: " << contig.multplty << std::endl;
        std::cout << std::endl; // Для разделения информации о каждом контиге
    }*/


    OverlapCollection overlap_collection = detect_adjacent_contigs(contig_collection, mink, maxk);
    
    /*for (const auto& pair : overlap_collection) {
        std::cout << "Key: " << pair.first << ", Value: ";
        for (const auto& overlap : pair.second) {
            std::cout << overlap.to_string() << " ";
        }
        std::cout << std::endl;
    }*/

    assign_multiplicity(contig_collection, overlap_collection);

    std::string outdpath = "output";
    write_summary(contig_collection, overlap_collection, filepath, outdpath);

    // Зафиксируем конечное время
    auto end_time = std::chrono::high_resolution_clock::now();

    // Вычислим время выполнения
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    // Выведем время выполнения в миллисекундах
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;


    return 0;
}
