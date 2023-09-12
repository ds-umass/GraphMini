//
// Created by ubuntu on 1/1/23.
//

#ifndef MINIGRAPH_CODEGEN_H
#define MINIGRAPH_CODEGEN_H

#include "common.h"
namespace minigraph
{

//    PlanIR create_plan(std::string pat, CodeGenConfig config, MetaData meta);
//    // PlanIR create_plan(std::string pat, CodeGenConfig config, MetaData meta);
//    std::string gen_code(PlanIR plan);

    /* brief Single Pattern Scheduling
     * pat: adjacency matrix of the input pattern
     * config: configuration of generated code
     * meta: metadata of the target graph
     * */
    std::string gen_code(const std::string& adj_mat, CodeGenConfig config, MetaData meta);
}

#endif //MINIGRAPH_CODEGEN_H