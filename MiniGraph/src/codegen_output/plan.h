#pragma once
#include "../backend/backend.h"
namespace minigraph
{
    using GraphType = Graph;
    using VertexSetType = VertexSet;
    void plan(const GraphType* graph, Context& ctx);
    uint64_t pattern_size();
}
