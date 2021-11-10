//
// Created by sunsh on 11/10/2021.
//

#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
    std::string path(argv[1]);
    std::ifstream file(path, std::ios::binary);

    uint32_t order_length;
    uint64_t embedding_count;
    file.read(reinterpret_cast<char *>(&order_length), sizeof(uint32_t));
    file.read(reinterpret_cast<char *>(&embedding_count), sizeof(uint64_t));

    uint32_t order_buffer_in_byte = order_length * sizeof(uint32_t);
    uint64_t embedding_buffer_in_byte = embedding_count * order_length * sizeof(uint32_t);

    uint32_t *matching_order = new uint32_t[order_length];
    uint32_t *embeddings = new uint32_t[embedding_count * order_length];

    file.read(reinterpret_cast<char *>(matching_order), order_buffer_in_byte);
    file.read(reinterpret_cast<char *>(embeddings), embedding_buffer_in_byte);


    /**
     * Output to terminal.
     */
    std::cout << "Matching Order Length: " << order_length << '\n';

    for (uint32_t i = 0; i < order_length; ++i) {
        std::cout << matching_order[i] <<' ';
    }
    std::cout <<'\n';

    std::cout << "Number of Embeddings: " << embedding_count << '\n';

    for (uint32_t i = 0; i < embedding_count; ++i) {
        uint32_t* embedding = embeddings + i * order_length;

        for (uint32_t j = 0; j < order_length; ++j) {
            std::cout << embedding[j] <<' ';
        }
        std::cout <<'\n';
    }

    delete[] matching_order;
    delete[] embeddings;

    return 0;
}