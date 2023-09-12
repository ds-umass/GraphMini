//
// Created by ubuntu on 2/14/23.
//
#include "logitem.h"
#include "logging.h"
#include <string>
#include <unordered_map>
#include <filesystem>
#include <fstream>
namespace minigraph
{
#define spliter ","
    void CompilerLog::save(std::string in_dir) {
        std::filesystem::path path = in_dir;
        path /= Constant::kExpCompileFile;
        bool file_exist = std::filesystem::is_regular_file(path);
        std::ofstream log(path, std::ios_base::app);
        if (!file_exist) {
            std::string items[] = {Constant::kExpID,
                                   Constant::kExpDataName,
                                   Constant::kExpPatternName,
                                   Constant::kExpPatternAdj, Constant::kExpPatternSize,
                                   Constant::kExpAdjMatType, Constant::kExpPruningType, Constant::kExpParallelType,
                                   Constant::kExpCodeGenTime, Constant::kExpCompilationTime,
                                   };
            // write header
            constexpr int num_items = std::size(items);
            for (int i = 0; i < num_items; i++) {
                log << items[i];
                if (i!=num_items-1) log << spliter;
            }
            log << "\n";
        }
        log << expId << spliter;
        log << dataName << spliter;
        log << patternName << spliter;
        log << patternAdj << spliter;
        log << patternSize << spliter;
        switch (adjMatType) {
            case AdjMatType::VertexInduced:
                log << "Vertex-Induced" << spliter;
                break;
            case AdjMatType::EdgeInduced:
                log << "Edge-Induced" << spliter;
                break;
            case AdjMatType::EdgeInducedIEP:
                log << "Edge-Induced-IEP" << spliter;
                break;
        }
        switch (pruningType) {
            case PruningType::None:
                log << "None" << spliter;
                break;
            case PruningType::Static:
                log << "Static" << spliter;
                break;
            case PruningType::Online:
                log << "Online" << spliter;
                break;
            case PruningType::CostModel:
                log << "CostModel" << spliter;
                break;
            case PruningType::Eager:
                log << "Eager" << spliter;
                break;
        }

        switch (parallelType) {
            case ParallelType::OpenMP:
                log << "OpenMP" << spliter;
                break;
            case ParallelType::TbbTop:
                log << "Tbb" << spliter;
                break;
            case ParallelType::Nested:
                log << "Nested" << spliter;
                break;
            case ParallelType::NestedRt:
                log << "NestedRt" << spliter;
                break;
            default:
                log << "Invalid ParallelType" << spliter;
                break;
        }
        log << codegenTime << spliter;
        log << compileTime;
        log << "\n";
        log.flush();
        log.close();
    }

    void RunnerLog::save(std::string in_dir) {
        std::filesystem::path path = in_dir;
        path /= Constant::kExpRunnerFile;
        bool file_exist = std::filesystem::is_regular_file(path);
        std::ofstream log(path, std::ios_base::app);
        if (!file_exist) {
            if (!file_exist) {
                std::string items[] = {Constant::kExpID,
                                       Constant::kExpNumThreads, Constant::kExpResult,
                                       Constant::kExpRunTime,
                                       Constant::kExpThreadMeanTime,Constant::kExpThreadMinTime,
                                       Constant::kExpThreadMaxTime, Constant::kExpThreadTimeSTD,
                                       Constant::kExpVertexAllocated, Constant::kExpVertexAllocatedPerThread,
                                       Constant::kExpMiniGraphAllocated, Constant::kExpMiniGraphAllocatedPerThread};
                // write header
                constexpr int num_items = std::size(items);
                for (int i = 0; i < num_items; i++) {
                    log << items[i];
                    if (i!=num_items-1) log << spliter;
                }
                log << "\n";
            }
        }

        log << expId << spliter;
        log << numThread << spliter;
        log << result << spliter;
        if (finished) log << runTime << spliter;
        else log << runTime << " (Time Out)" << spliter;
        log << threadMeanTime << spliter;
        log << threadMinTime << spliter;
        log << threadMaxTime << spliter;
        log << threadTimeSTD << spliter;
        log << ToReadableSize(vertexAllocated) << spliter;
        log << ToReadableSize(vertexAllocated / numThread) << spliter;
        log << ToReadableSize(miniGraphAllocated) << spliter;
        log << ToReadableSize(miniGraphAllocated / numThread);
        log << "\n";
        log.flush();
        log.close();
    }
}