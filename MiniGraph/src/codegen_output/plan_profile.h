//
// Created by ubuntu on 2/24/23.
//

#ifndef MINIGRAPH_PLAN_PROFILE_H
#define MINIGRAPH_PLAN_PROFILE_H
#include "../backend_prof/backend.h"
namespace minigraph
{
    using GraphType = Graph;
    using VertexSetType = VertexSet;
    void plan(const GraphType* graph, Context& ctx);
    uint64_t pattern_size();
}

#endif //MINIGRAPH_PLAN_PROFILE_H
