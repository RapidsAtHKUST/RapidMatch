#ifndef SUBGRAPHMATCHING_ENCODER_H
#define SUBGRAPHMATCHING_ENCODER_H

#include <relation/catalog.h>

enum RelationStructure {
    EncodedTrieRelation,
    HashRelation,
    TrieRelation
};

class encoder {
public:
    double encoding_time_;
    double build_relation_time_;
    double projection_time_;
    double convert_to_sparse_bp_time_;

private:
    uint32_t max_vertex_id_;
    uint32_t* temp_buffer1_;
    uint32_t* temp_buffer2_;
    const Graph* query_graph_;
    uint32_t* join_plan_;

private:
    void convert_to_trie_relation(catalog *storage);
    void convert_to_trie_relation(catalog *storage, uint32_t u, uint32_t v);

    void convert_to_hash_relation(catalog *storage);
    void convert_to_hash_relation(catalog *storage, uint32_t u, uint32_t v);

    void convert_to_encoded_relation(catalog *storage);
    void convert_to_encoded_relation(catalog *storage, uint32_t u, uint32_t v);

    void convert_trie_relation_to_sparse_bitmap(catalog *storage);
    void convert_encoded_relation_to_sparse_bitmap(catalog *storage);
    void convert_hash_relation_to_sparse_bitmap(catalog *storage);

public:
    encoder(const Graph* query_graph, uint32_t max_vertex_id) :
            encoding_time_(0), build_relation_time_(0), projection_time_(0), convert_to_sparse_bp_time_(0),
            max_vertex_id_(max_vertex_id) {
        query_graph_ = query_graph;
        temp_buffer1_ = new uint32_t[max_vertex_id_];
        temp_buffer2_ = new uint32_t[max_vertex_id_];
        memset(temp_buffer1_, 0, max_vertex_id_ * sizeof(uint32_t));
    }

    ~encoder() {
        delete[] temp_buffer1_;
        delete[] temp_buffer2_;
    }

    void execute(catalog *storage, RelationStructure relation_type, bool enable_sparsebp, uint32_t *join_plan);

    void print_metrics();
};


#endif //SUBGRAPHMATCHING_ENCODER_H
