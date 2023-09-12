//
// Created by ubuntu on 2/14/23.
//

#ifndef MINIGRAPH_LOGITEM_H
#define MINIGRAPH_LOGITEM_H
#include "constant.h"
#include "typedef.h"
#include <string>
namespace minigraph {
    struct CompilerLog
    {
        int expId{0};
        int patternSize{0};
        double compileTime{0.0};
        double codegenTime{0.0};
        PruningType pruningType;
        ParallelType parallelType;
        AdjMatType adjMatType;
        std::string patternName;
        std::string patternAdj;
        std::string dataName;
        void save(std::string in_dir);
//        void read(std::string in_dir, int exp_id);
    };

    struct RunnerLog {
        int expId{0};
        bool finished{true};
        long long result;
        long long numThread;
        long long vertexAllocated;
        long long miniGraphAllocated;
        double runTime;
        double threadMeanTime;
        double threadMinTime;
        double threadMaxTime;
        double threadTimeSTD;
        void save(std::string in_dir);
//        void read(std::string in_dir, int exp_id);
    };
}
#endif //MINIGRAPH_LOGITEM_H
