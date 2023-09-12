////
//// Created by ubuntu on 1/18/23.
////
//
//#include "constant.h"
//#include <string>
//#include <limits>
//
//namespace minigraph
//{
//    const std::string Constant::kMetaFile = "meta.txt";
//    const std::string Constant::kMetaNumVertex = "NUM_VERTEX";
//    const std::string Constant::kMetaNumEdge = "NUM_EDGE";
//    const std::string Constant::kMetaNumTriangle = "NUM_TRIANGLE";
//    const std::string Constant::kMetaMaxDegree = "MAX_DEGREE";
//    const std::string Constant::kMetaMaxOffset = "MAX_OFFSET";
//    const std::string Constant::kMetaMaxTriangle = "MAX_TRIANGLE";
//    // Graph Data
//    std::string Constant::kDataFile = "snap.txt";
//    const std::string Constant::kIndptrU64File = "indptr_u64.bin";
//    const std::string Constant::kIndicesU64File = "indices_u64.bin";
//    const std::string Constant::kTriangleU64File = "triangle_u64.bin";
//    const std::string Constant::kOffsetU64File = "offset_u64.bin";
//    const std::string Constant::kDegreeU64File = "degree_u64.bin";
//
//    const std::string Constant::kIndptrU32File = "indptr_u32.bin";
//    const std::string Constant::kIndicesU32File = "indices_u32.bin";
//    const std::string Constant::kTriangleU32File = "triangle_u32.bin";
//    const std::string Constant::kOffsetU32File = "offset_u32.bin";
//    const std::string Constant::kDegreeU32File = "degree_u32.bin";
//
//    // Profiling Files
//    const std::string Constant::kExpCompileFile = "compile_log.csv";
//    const std::string Constant::kExpRunnerFile = "run_log.csv";
//    const std::string Constant::kProfileFile = "profile.csv";
//    const std::string Constant::kProfileMetaFile = "profile_meta.txt";
//    const std::string Constant::kProfileMetaNumVertex = "NUM_VERTEX";
//    // Experiment Configuration
//    const std::string Constant::kExpPatternName = "PATTERN_NAME";
//    const std::string Constant::kExpID = "EXP_ID";
//    const std::string Constant::kExpDataName = "DATA_NAME";
//     const std::string Constant::kExpPatternAdj = "PATTERN_ADJ";
//     const std::string Constant::kExpPatternSize = "PATTERN_SIZE";
//     const std::string Constant::kExpAdjMatType = "ADJMAT_TYPE";
//     const std::string Constant::kExpPruningType = "PRUNING_TYPE";
//    // Experiment Items -- Compiler
//    const std::string Constant::kExpParallelType = "PARALLEL_TYPE";
//    const std::string Constant::kExpCodeGenTime = "CODEGEN_TIME";
//    const std::string Constant::kExpCompilationTime = "COMPILE_TIME";
//
//    // Experiment Items -- Runner
//    const std::string Constant::kExpRunTime = "RUN_TIME";
//    const std::string Constant::kExpResult = "RESULT";
//    const std::string Constant::kExpNumThreads = "NUMBER_THREAD";
//    const std::string Constant::kExpThreadMeanTime = "THREAD_MEAN_TIME";
//    const std::string Constant::kExpThreadMinTime = "THREAD_MIN_TIME";
//    const std::string Constant::kExpThreadMaxTime = "THREAD_MAX_TIME";
//    const std::string Constant::kExpThreadTimeSTD = "THREAD_TIME_STD";
//
//    const std::string Constant::kExpVertexAllocated = "VERTEX_ALLOCATED";
//    const std::string Constant::kExpMiniGraphAllocated = "MINIGRAPH_ALLOCATED";
//    const std::string Constant::kExpVertexAllocatedPerThread = "VERTEX_ALLOCATED_PER_THREAD";
//    const std::string Constant::kExpMiniGraphAllocatedPerThread = "MINIGRAPH_ALLOCATED_PER_THREAD";
//
//    const std::string Constant::kSetCompPerLoopPerVertexBin = "VID_TO_SET_COMP.bin";
//    const std::string Constant::kGmCompPerLoopPerVertexBin = "gm_comp_perloop_pervertex.bin";
//    const std::string Constant::kAdjReadPerLoopPerVertexBin = "VID_TO_FREQ.bin";
//    const std::string Constant::kSetCompPerLoopPerVertexCsv = "VID_TO_SET_COMP.csv";
//    const std::string Constant::kGmCompPerLoopPerVertexCsv = "gm_comp_perloop_pervertex.csv";
//    const std::string Constant::kAdjReadPerLoopPerVertexCsv = "VID_TO_FREQ.csv";
//}