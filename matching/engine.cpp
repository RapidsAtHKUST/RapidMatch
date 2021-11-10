#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <computesetintersection.h>

#include "matchingcommand.h"
#include "graph/graph.h"
#include "relation/catalog.h"
#include "execution_tree_generator.h"
#include "pretty_print.h"
#include "global_variables.h"
#include "preprocessor.h"
#include "encoder.h"
#include "query_plan_generator.h"

void execute_within_time_limit(execution_tree* tree, catalog* catalog, uint64_t output_limit, uint64_t time_limit) {
    g_exit = false;
    std::future<uint64_t> future = std::async(std::launch::async, [tree, catalog, output_limit](){
        return tree->execute(*catalog, output_limit);;
    });

    std::future_status status;
    do {
        status = future.wait_for(std::chrono::seconds(time_limit));
        if (status == std::future_status::deferred) {
            std::cout << "Deferred\n";
            exit(-1);
        } else if (status == std::future_status::timeout) {
            g_exit = true;
        }
    } while (status != std::future_status::ready);
}

void output_embedding(const std::string& path, uint32_t *matching_order, uint32_t order_length,
                      uint32_t *embeddings, uint64_t embedding_count) {
    std::cout << "Embedding Number Limit: " << OUTPUT_RESULT_NUM_LIMIT << std::endl;
    std::cout << "Number of Embeddings Dumped to File: " << embedding_count << std::endl;

    std::ofstream embedding_file(path, std::ios::binary);

    uint32_t order_size = order_length * sizeof(uint32_t);
    uint64_t embedding_size = embedding_count * order_length * sizeof(uint32_t);
    // Output the size of the matching order.
    embedding_file.write(reinterpret_cast<const char *>(&order_length), sizeof(uint32_t));
    // output the size of the embeddings.
    embedding_file.write(reinterpret_cast<const char *>(&embedding_count), sizeof(uint64_t));

    embedding_file.write(reinterpret_cast<const char *>(matching_order), order_size);

    embedding_file.write(reinterpret_cast<const char *>(embeddings), embedding_size);
}

int main(int argc, char** argv) {
    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string input_order_type = command.getOrderType();
    std::string input_export_plan_path = command.getExportPlanPath();
    std::string input_max_embedding_num = command.getMaximumEmbeddingNum();
    std::string input_time_limit = command.getTimeLimit();
    std::string input_order_num = command.getOrderNum();
    std::string input_csr_file_path = command.getCSRFilePath();
    std::string input_import_plan_path = command.getImportPlanPath();
    std::string input_order = command.getInputOrder();
    std::string input_enable_preprocessor = command.getPreprocessor();
    std::string input_output_path = command.getOutputPath();
    /**
     * Output the command line information.
     */
    std::cout << "Command Line:" << std::endl;
    std::cout << "\tData Graph CSR: " << input_csr_file_path << std::endl;
    std::cout << "\tData Graph: " << input_data_graph_file << std::endl;
    std::cout << "\tQuery Graph: " << input_query_graph_file << std::endl;
    std::cout << "\tOrder Type: " << input_order_type << std::endl;
    std::cout << "\tExport Plan Path: " << input_export_plan_path << std::endl;
    std::cout << "\tImport Plan Path: " << input_import_plan_path << std::endl;
    std::cout << "\tOutput Limit: " << input_max_embedding_num << std::endl;
    std::cout << "\tTime Limit (seconds): " << input_time_limit << std::endl;
    std::cout << "\tOrder Num: " << input_order_num << std::endl;
    std::cout << "\tInput Order: " << input_order << std::endl;
    std::cout << "\tEnable Preprocessor: " << input_enable_preprocessor << std::endl;
#ifdef ENABLE_OUTPUT
    std::cout << "\tOutput Path: " << input_output_path << std::endl;
#endif
    std::cout << "--------------------------------------------------------------------" << std::endl;

    /**
     * Output the configuration.
     */
    std::cout << "Configuration:" << std::endl;
#if RELATION_STRUCTURE == 0
    RelationStructure relation_type = RelationStructure::EncodedTrieRelation;
    std::cout << "\tRelation Structure: Encoded trie relation" << std::endl;
#elif RELATION_STRUCTURE == 1
    RelationStructure relation_type = RelationStructure::HashRelation;
    std::cout << "\tRelation Structure: Hash relation" << std::endl;
#elif RELATION_STRUCTURE == 2
    RelationStructure relation_type = RelationStructure::TrieRelation;
    std::cout << "\tRelation Structure: Trie relation" << std::endl;
#endif

#ifdef HOMOMORPHISM
    std::cout << "\tEmbedding Structure: Subgraph Homomorphism " << std::endl;
#else
    std::cout << "\tEmbedding Structure: Subgraph Isomorphism " << std::endl;
#endif

#ifdef COUNT_RESULTS
    std::cout << "\tCount Results: True" << std::endl;
#else
    std::cout << "\tCount Results: False" << std::endl;
#endif

#if HYBRID == 0
    std::cout << "\tHybrid Set Intersection: Enable" << std::endl;
#else
    std::cout << "\tHybrid Set Intersection: Disable" << std::endl;
#endif

#if SI == 0
    std::cout << "\tSIMD Set Intersection: AVX2" << std::endl;
#elif SI == 1
    std::cout << "\tSIMD Set Intersection: AVX512" << std::endl;
#elif SI == 2
    std::cout << "\tSIMD Set Intersection: NONE" << std::endl;
#endif

    std::cout << "\tHash Table Ratio of Pair-wise Join: " << HASH_TABLE_RATIO << std::endl;

#ifdef INTERSECTION_CACHE
    std::cout << "\tIntersection Cache: Enabled" << std::endl;
#else
    std::cout << "\tIntersection Cache: Disabled" << std::endl;
#endif

#ifdef FAILING_SET_PRUNING
    std::cout << "\tFailing Set Pruning: Enabled" << std::endl;
#else
    std::cout << "\tFailing Set Pruning: Disabled" << std::endl;
#endif

#ifdef OUTPUT_OPTIMIZATION
    std::cout << "\tOutput Optimization: Enabled" << std::endl;
#else
    std::cout << "\tOutput Optimization: Disabled" << std::endl;
#endif

#ifdef SPARSE_BITMAP
    bool enable_sparsebp = true;
    std::cout << "\tSparse Bitmap: Enabled" << std::endl;
#else
    bool enable_sparsebp = false;
    std::cout << "\tSparse Bitmap: Disabled" << std::endl;
#endif

    std::cout << "--------------------------------------------------------------------" << std::endl;
    /**
     * Load input graphs.
     */
    std::cout << "Load graphs..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->buildCoreTable();

    auto data_graph = new Graph(true);

    if (input_csr_file_path.empty()) {
        data_graph->loadGraphFromFile(input_data_graph_file);
    }
    else {
        std::string degree_file_path = input_csr_file_path + "_deg.bin";
        std::string edge_file_path = input_csr_file_path + "_adj.bin";
        std::string label_file_path = input_csr_file_path + "_label.bin";
        data_graph->loadGraphFromFileCompressed(degree_file_path, edge_file_path, label_file_path);
    }

    auto end = std::chrono::high_resolution_clock::now();

    double load_graphs_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

#ifdef ENABLE_OUTPUT
    output_buffer = new uint32_t[OUTPUT_RESULT_NUM_LIMIT * (uint64_t)query_graph->getVerticesCount()];
#endif
    std::cout << "-----" << std::endl;
    std::cout << "Query Graph Meta Information" << std::endl;
    query_graph->printGraphMetaData();
    std::cout << "-----" << std::endl;
    std::cout << "Data Graph Meta Information" << std::endl;
    data_graph->printGraphMetaData();

    std::cout << "--------------------------------------------------------------------" << std::endl;

    std::cout << "Preprocess..." << std::endl;

    double query_time = 0;
    bool enable_preprocessor = true;
    if (input_enable_preprocessor == "false") {
        enable_preprocessor = false;
    }
    else if (input_enable_preprocessor == "true") {
        enable_preprocessor = true;
    }
    else {
        std::cout << "The argument of -preprocess is not supported." << std::endl;
        exit(-1);
    }

    // Execute preprocessor.
    auto pp = new preprocessor();
    auto storage = new catalog(query_graph, data_graph);

    pp->execute(query_graph, data_graph, storage, enable_preprocessor);
    pp->print_metrics();
    query_time += pp->preprocess_time_;

    uint32_t order_num = 0;
    if (input_order_num == "MAX") {
        order_num = std::numeric_limits<uint32_t>::max();
    }
    else {
        sscanf(input_order_num.c_str(), "%u", &order_num);
    }


    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Generate query plan..." << std::endl;

    std::vector<std::vector<uint32_t>> spectrum;
    if (input_order_type == "Spectrum") {
        execution_tree_generator::generate_wcoj_orders(query_graph, spectrum, order_num);
    }
    else if (input_order_type == "Import") {
    }
    else if (input_order_type == "nd") {
        query_plan_generator::generate_query_plan_with_nd(query_graph, storage, spectrum);
        query_plan_generator::print_metrics();
    }
    else if (input_order_type == "Test") {
        std::vector<uint32_t> order;
        query_plan_generator::generate_query_plan_for_test(query_graph, order);
        spectrum.emplace_back(order);
        std::vector<std::vector<uint32_t>> temp;
        query_plan_generator::print_vertex_orders(query_graph, spectrum, temp);
    }
    else if (input_order_type == "Input") {
        std::vector<uint32_t> order;
        std::stringstream ss(input_order);

        for (uint32_t u; ss >> u;) {
            order.push_back(u);

            if (ss.peek() == ',')
                ss.ignore();
        }

        spectrum.emplace_back(order);
        std::vector<std::vector<uint32_t>> temp;
        query_plan_generator::print_vertex_orders(query_graph, spectrum, temp);
    }
    else {
        std::cout << "The specified order type '" << input_order_type << "' is not supported." << std::endl;
        exit(-1);
    }
    query_time += query_plan_generator::ordering_time_;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Encode..." << std::endl;
    auto en = new encoder(query_graph, data_graph->getVerticesCount());
    en->execute(storage, relation_type, enable_sparsebp, spectrum.back().data());
    en->print_metrics();
    query_time += en->encoding_time_;

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Catalog info..." << std::endl;
    storage->print_metrics();

    std::cout << "--------------------------------------------------------------------" << std::endl;

    std::cout << "Enumerate..." << std::endl;
    uint64_t output_limit = 0;
    uint64_t embedding_count = 0;
    if (input_max_embedding_num == "MAX") {
        output_limit = std::numeric_limits<uint64_t>::max();
    }
    else {
        sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
    }

    double enumeration_time_in_ns = 0;
    uint64_t call_count = 0;
    uint64_t time_limit = 0;
    if (input_time_limit == "MAX") {
        time_limit = std::numeric_limits<uint64_t>::max();
    }
    else {
        sscanf(input_time_limit.c_str(), "%zu", &time_limit);
    }

    if (input_order_type != "Import") {
        std::cout << "Number of orders: " << spectrum.size() << std::endl;

        bool output_query_plan = false;
        std::ofstream ofs;
        if (!input_export_plan_path.empty()) {
            output_query_plan = true;
            ofs.open(input_export_plan_path);
        }
        uint32_t order_id = 0;
        for (auto &order : spectrum) {
            std::cout << "------------------------" << std::endl;
            std::cout << "Spectrum Order " << ++order_id << " Matching order: "<< order << std::endl;

            auto tree = execution_tree_generator::generate_single_node_execution_tree(order);
            if (output_query_plan) {
                ofs << tree->export_plan() << std::endl;
            }

            execute_within_time_limit(tree, storage, output_limit, time_limit);

            embedding_count = tree->get_output_count();
            enumeration_time_in_ns = tree->get_execution_time();

            printf("Spectrum Order %u Enumerate time (seconds): %.6lf\n", order_id,
                   NANOSECTOSEC(enumeration_time_in_ns));
            printf("Spectrum Order %u #Embeddings: %zu\n", order_id, embedding_count);
            printf("Spectrum Order %u Throughput per second: %.6lf\n", order_id, embedding_count / NANOSECTOSEC(enumeration_time_in_ns));
            printf("Spectrum Order %u Call Count: %zu\n", order_id, call_count);
            printf("Spectrum Order %u (vertex, pr_count, si_count, si_cost): %s\n", order_id, tree->get_cost_metrics().c_str());
            printf("Spectrum Order %u Invalid Partial Results Count (total, core, leaf): %s\n",  order_id, tree->get_invalid_pr_count().c_str());
            printf("Spectrum Order %u Fail Count (total, si_fail, iso_conflict): %s\n",  order_id, tree->get_fail_count().c_str());
            printf("Spectrum Order %u Merge, Galloping: %zu, %zu\n", order_id, ComputeSetIntersection::merge_cnt_, ComputeSetIntersection::galloping_cnt_);
            ComputeSetIntersection::galloping_cnt_ = 0;
            ComputeSetIntersection::merge_cnt_ = 0;
            delete tree;
        }
    }
    else {
        std::ifstream ifs(input_import_plan_path);
        std::string query_pan((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
        auto* tree = new execution_tree(query_pan);
        std::cout << "Query plan:\n" << tree->export_plan() << std::endl;
        execute_within_time_limit(tree, storage, output_limit, time_limit);
        embedding_count = tree->get_output_count();
        enumeration_time_in_ns = tree->get_execution_time();

        delete tree;
    }

#ifdef ENABLE_OUTPUT
    std::cout << "Dump results to file..." << std::endl;
    output_embedding(input_output_path, spectrum[0].data(), spectrum[0].size(),
                     output_buffer, embedding_count);
    delete[] output_buffer;
#endif

    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;
    /**
     * Release the allocated memories.
     */

    delete query_graph;
    delete data_graph;
    delete storage;
    delete en;
    delete pp;

    /**
     * End.
     */
    std::cout << "--------------------------------------------------------------------" << std::endl;
    query_time += enumeration_time_in_ns;
    printf("Load graphs time (seconds): %.6lf\n", NANOSECTOSEC(load_graphs_time_in_ns));
    printf("Query time (seconds): %.6lf\n", NANOSECTOSEC(query_time));
    printf("#Embeddings: %zu\n", embedding_count);

    std::cout << "End." << std::endl;
    return 0;
}
