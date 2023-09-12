//
// Created by ubuntu on 1/19/23.
//

#ifndef MINIGRAPH_UTILITY_H
#define MINIGRAPH_UTILITY_H
#include "graph.h"
#include "vertex_set.h"
#include <iosfwd>
namespace minigraph
{
    template<typename IdType>
    Graph<IdType>* load_bin(std::string in_dir, bool _mmap);

    template<typename IdType>
    std::ostream &operator<<(std::ostream &os, const VertexSet<IdType> &dt);
}
#endif //MINIGRAPH_UTILITY_H
