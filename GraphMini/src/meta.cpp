//
// Created by ubuntu on 1/18/23.
//

#include "meta.h"
#include "constant.h"
#include <assert.h>
#include <iosfwd>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <iterator>
#include <sstream>
namespace minigraph
{

    void MetaData::save(std::string ind_dir) {
        std::filesystem::path path = ind_dir;
        path /= Constant::kMetaFile;
        std::ofstream meta(path, std::ios_base::out);
        meta << Constant::kMetaNumVertex << "\t" << num_vertex << "\n";
        meta << Constant::kMetaNumEdge << "\t" << num_edge << "\n";
        meta << Constant::kMetaNumTriangle << "\t" << num_triangle << "\n";
        meta << Constant::kMetaMaxDegree << "\t" << max_degree << "\n";
        meta << Constant::kMetaMaxOffset << "\t" << max_offset << "\n";
        meta << Constant::kMetaMaxTriangle << "\t" << max_triangle << "\n";
        meta.close();
    }

    void MetaData::read(std::string in_dir) {
        std::filesystem::path path = in_dir;
        path /= Constant::kMetaFile;
        assert(std::filesystem::is_regular_file(path) && "Meta file does not exists" );
        std::ifstream meta(path);
        std::string line;
        std::unordered_map<std::string, uint64_t> items;
        while (std::getline(meta, line)) {
            std::istringstream iss(line);
            std::vector<std::string> kv{std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}};
            if (kv.size() < 2) {
                break;
            }
            items[kv[0]] = std::stoull(kv[1]);
        }
        assert(items.count(Constant::kMetaNumVertex) > 0);
        assert(items.count(Constant::kMetaNumEdge) > 0);
        assert(items.count(Constant::kMetaNumTriangle) > 0);
        assert(items.count(Constant::kMetaMaxDegree) > 0);
        assert(items.count(Constant::kMetaMaxOffset) > 0);
        assert(items.count(Constant::kMetaMaxTriangle) > 0);
        num_vertex = items[Constant::kMetaNumVertex];
        num_edge = items[Constant::kMetaNumEdge];
        num_triangle = items[Constant::kMetaNumTriangle];
        max_degree = items[Constant::kMetaMaxDegree];
        max_offset = items[Constant::kMetaMaxOffset];
        max_triangle = items[Constant::kMetaMaxTriangle];
        num_triangle /= 6; // remove automorphism to make it in consistent with GraphPi
    }
}