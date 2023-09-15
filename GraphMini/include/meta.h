//
// Created by ubuntu on 1/18/23.
//

#ifndef MINIGRAPH_META_H
#define MINIGRAPH_META_H
#include <string>
#include <stdint.h>

namespace minigraph
{
    class MetaData {
    public:
        uint64_t num_vertex{0};
        uint64_t num_edge{0};
        uint64_t num_triangle{0};
        uint64_t max_degree{0};
        uint64_t max_offset{0};
        uint64_t max_triangle{0};

        MetaData() = default;
        MetaData(uint64_t _num_vertex, uint64_t _num_edge, uint64_t _num_triangle,
                 uint64_t _max_degree, uint64_t _max_offset, uint64_t _max_triangle):
                num_vertex{_num_vertex}, num_edge{_num_edge}, num_triangle{_num_triangle},
                max_degree{_max_degree}, max_offset{_max_offset}, max_triangle{_max_triangle}{};

        void save(std::string in_dir);
        void read(std::string in_dir);;
    };
}
#endif //MINIGRAPH_META_H
