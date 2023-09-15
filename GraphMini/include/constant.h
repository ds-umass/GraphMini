//
// Created by ubuntu on 1/18/23.
//

#ifndef MINIGRAPH_CONSTANT_H
#define MINIGRAPH_CONSTANT_H
#include <string>
#include <limits>
namespace minigraph
{
    class Constant {
    public:
        // Constant
        static constexpr size_t kResizeIncrement = 4096;
        static constexpr double kAllocScale = 1.25;
        static constexpr size_t kAllocNoScale = 1;
        static constexpr size_t kBufferSize = 64;
        static constexpr size_t kGigabytes = 1 * 1024 * 1024 * 1024;
        static constexpr size_t kMegabytes = 1 * 1024 * 1024;
        static constexpr size_t kKilobytes = 1 * 1024;

        template<typename T>
        static constexpr T EmptyID() {
            return std::numeric_limits<T>::max();
        };

        // Graph Meta Data
        inline static const std::string kMetaFile = "meta.txt";
        inline static const std::string kMetaNumVertex = "NUM_VERTEX";
        inline static const std::string kMetaNumEdge = "NUM_EDGE";
        inline static const std::string kMetaNumTriangle = "NUM_TRIANGLE";
        inline static const std::string kMetaMaxDegree = "MAX_DEGREE";
        inline static const std::string kMetaMaxOffset = "MAX_OFFSET";
        inline static const std::string kMetaMaxTriangle = "MAX_TRIANGLE";
        // Graph Data
        inline static const std::string kDataFile = "snap.txt";
        inline static const std::string kIndptrU64File = "indptr_u64.bin";
        inline static const std::string kIndicesU64File = "indices_u64.bin";
        inline static const std::string kTriangleU64File = "triangle_u64.bin";
        inline static const std::string kOffsetU64File = "offset_u64.bin";
        inline static const std::string kDegreeU64File = "degree_u64.bin";

        inline static const std::string kIndptrU32File = "indptr_u32.bin";
        inline static const std::string kIndicesU32File = "indices_u32.bin";
        inline static const std::string kTriangleU32File = "triangle_u32.bin";
        inline static const std::string kOffsetU32File = "offset_u32.bin";
        inline static const std::string kDegreeU32File = "degree_u32.bin";

        // Profiling Files
        inline static const std::string kExpCompileFile = "compile_log.csv";
        inline static const std::string kExpRunnerFile = "run_log.csv";
        inline static const std::string kProfileFile = "compile_and_run.csv";
        inline static const std::string kProfileMetaFile = "profile_meta.txt";
        inline static const std::string kProfileMetaNumVertex = "NUM_VERTEX";
        // Experiment Configuration
        inline static const std::string kExpPatternName = "PATTERN_NAME";
        inline static const std::string kExpID = "EXP_ID";
        inline static const std::string kExpDataName = "DATA_NAME";
        inline static const std::string kExpPatternAdj = "PATTERN_ADJ";
        inline static const std::string kExpPatternSize = "PATTERN_SIZE";
        inline static const std::string kExpAdjMatType = "ADJMAT_TYPE";
        inline static const std::string kExpPruningType = "PRUNING_TYPE";
        // Experiment Items -- Compiler
        inline static const std::string kExpParallelType = "PARALLEL_TYPE";
        inline static const std::string kExpCodeGenTime = "CODEGEN_TIME";
        inline static const std::string kExpCompilationTime = "COMPILE_TIME";

        // Experiment Items -- Runner
        inline static const std::string kExpRunTime = "RUN_TIME";
        inline static const std::string kExpResult = "RESULT";
        inline static const std::string kExpNumThreads = "NUMBER_THREAD";
        inline static const std::string kExpThreadMeanTime = "THREAD_MEAN_TIME";
        inline static const std::string kExpThreadMinTime = "THREAD_MIN_TIME";
        inline static const std::string kExpThreadMaxTime = "THREAD_MAX_TIME";
        inline static const std::string kExpThreadTimeSTD = "THREAD_TIME_STD";

        inline static const std::string kExpVertexAllocated = "VERTEX_ALLOCATED";
        inline static const std::string kExpMiniGraphAllocated = "MINIGRAPH_ALLOCATED";
        inline static const std::string kExpVertexAllocatedPerThread = "VERTEX_ALLOCATED_PER_THREAD";
        inline static const std::string kExpMiniGraphAllocatedPerThread = "MINIGRAPH_ALLOCATED_PER_THREAD";

        inline static const std::string kSetCompPerLoopPerVertexBin = "VID_TO_SET_COMP.bin";
        inline static const std::string kGmCompPerLoopPerVertexBin = "gm_comp_perloop_pervertex.bin";
        inline static const std::string kAdjReadPerLoopPerVertexBin = "VID_TO_FREQ.bin";
        inline static const std::string kSetCompPerLoopPerVertexCsv = "VID_TO_SET_COMP.csv";
        inline static const std::string kGmCompPerLoopPerVertexCsv = "gm_comp_perloop_pervertex.csv";
        inline static const std::string kAdjReadPerLoopPerVertexCsv = "VID_TO_FREQ.csv";
    };
}
#endif //MINIGRAPH_CONSTANT_H
