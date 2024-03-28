#include <iostream>
#include <string>

#include "contigs.hpp"
#include "overlaps.hpp"
#include "assign_multiplicity.hpp"


int main() {
    std::string filepath = "/home/nazar/projects/CombinatorCPP/test_contigs_a5_0.fasta";
    int maxk=25;
    int mink=2;
    
    // Получение коллекции контигов
    ContigCollection contig_collection = get_contig_collection(filepath, maxk);

    OverlapCollection overlap_collection = detect_adjacent_contigs(contig_collection, mink, maxk);

    // Output the results
    std::cout << "Overlap Collection:" << std::endl;
    std::cout << overlap_collection.to_string() << std::endl;

    assign_multiplicity(contig_collection, overlap_collection);

    for (const auto& contig : contig_collection) {
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
    }


    return 0;
}