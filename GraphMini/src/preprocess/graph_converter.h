//
// Created by ubuntu on 1/1/23.
//

#ifndef MINIGRAPH_GRAPH_CONVERTER_H
#define MINIGRAPH_GRAPH_CONVERTER_H
#include "typedef.h"
#include <stdint.h>
#include <filesystem>
#include <vector>
namespace minigraph {
    class GraphConverter {
    private:
        std::vector<uint64_t> indices, indptr, triangles, offsets, degrees;
        uint64_t v_num{0}, e_num{0}, tri_num{0}, max_deg{0}, max_offset{0}, max_tri{0};
        std::filesystem::path in_dir;
        std::filesystem::path data_file;
        void load_txt();
        void save_meta();
        // save indices to unsigned 32-bit integer format
        void save_bin_u32();

        // save indices to unsigned 64-bit integer format
        void save_bin();

    public:
        GraphConverter() = default;

        void convert(std::filesystem::path input_dir);

    };
}
#endif //MINIGRAPH_GRAPH_CONVERTER_H
