//
// Created by ubuntu on 1/2/23.
//

#ifndef MINIGRAPH_TYPEDEF_H
#define MINIGRAPH_TYPEDEF_H

namespace minigraph {
//---------- Basic Types ----------
    enum AdjMatType {
       VertexInduced = 0, EdgeInduced = 1, EdgeInducedIEP = 2,
    };

    enum class PruningType {
        None = 0, // not using minigraph
        Static = 1, // use iff not introducing any redundant set operation
        Eager = 2, // use for all potential cases
        Online = 3, // use iff pruned adj will be used once
        CostModel = 4, // lazy + cost model
    };

    enum class ParallelType {
        OpenMP = 0,
        TbbTop = 1, // only parallel at top level
        Nested = 2, // aggressively nested
        NestedRt = 3, // decide whether to parallel at codegen_output
//        Distributed =4 // TODO: Implement it
    };

    enum class RunnerType {
        Benchmark,
        Profiling // TODO: Implement it
    };

    struct CodeGenConfig {
        AdjMatType adjMatType = AdjMatType::VertexInduced;
        PruningType pruningType = PruningType::None;
        ParallelType parType = ParallelType::OpenMP;
        RunnerType runnerType = RunnerType::Benchmark;
    };


    // ---------- Utility Functions & Structs ----------
    struct Counter {
        struct value_type {
            template<typename T>
            value_type(const T &) {}
        };
        void push_back(const value_type &) { ++count; }
        unsigned long count {0};
    };
}
#endif //MINIGRAPH_TYPEDEF_H
