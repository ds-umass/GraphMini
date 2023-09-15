//
// Created by ubuntu on 1/1/23.
//
#include "codegen.h"
#include "logging.h"
#include "configure.h"
#include "common.h"
#include <vector>
#include <iostream>
#include <fmt/format.h>
#include <filesystem>
#include <fstream>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <array>
#include <chrono>
#include <thread>

std::string exec(const char *cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        return "runtime error\n";
        // throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

using namespace minigraph;

struct AppConfig {
    int exp_id;
    CodeGenConfig codegen;
    std::string data_name;
    std::string pattern_name;
    std::string pat;
};

std::filesystem::path output_dir(AppConfig config) {
    std::filesystem::path out = PROJECT_LOG_DIR;
    out /= config.data_name;
    out /= config.pattern_name;
    if (config.codegen.adjMatType == AdjMatType::EdgeInduced) {
        out /= "edge_induced";
    } else if (config.codegen.adjMatType == AdjMatType::VertexInduced) {
        out /= "vertex_induced";
    } else {
        out /= "edge_induced_iep";
    };

    return out;
}

std::filesystem::path output_log(AppConfig config) {
    auto out = output_dir(config);
    std::string prefix;
    if (config.codegen.parType == ParallelType::OpenMP) prefix = "omp_";
    if (config.codegen.parType == ParallelType::TbbTop) prefix = "tbb_top_";
    if (config.codegen.parType == ParallelType::Nested) prefix = "tbb_nested_";
    if (config.codegen.parType == ParallelType::NestedRt) prefix = "tbb_nested_rt_";
    if (config.codegen.pruningType == PruningType::None) {
        out /= prefix + "baseline.txt";
    } else if (config.codegen.pruningType == PruningType::Eager) {
        out /= prefix + "eager.txt";
    } else if (config.codegen.pruningType == PruningType::Static) {
        out /= prefix + "lazy.txt";
    } else if (config.codegen.pruningType == PruningType::Online) {
        out /= prefix + "online.txt";
    }
    return out;
}

std::filesystem::path code_path() {
    std::filesystem::path code_file(PROJECT_SOURCE_DIR);
    code_file /= "src";
    code_file /= "codegen_output";
    code_file /= "plan_profile.cpp";
    return code_file;
}

std::filesystem::path get_data_dir(AppConfig config) {
    std::filesystem::path dd = "/home/ubuntu/";
    dd /= "dataset";
    dd /= config.data_name;
    return dd;
}

void compile(AppConfig config) {
    CompilerLog log;
    LOG(INFO) << "START EXPERIMENT";
    LOG(INFO) << "Pattern Name\t" << config.pattern_name;
    LOG(INFO) << "AdjMat\t" << config.pat;
    LOG(INFO) << "Graph\t" << config.data_name;
    MetaData meta;
    meta.read(get_data_dir(config));
    Timer t;
    std::string code = gen_code(config.pat, config.codegen, meta);
    auto codegen_t = t.Passed();
    t.Reset();
    std::ofstream out_file(code_path());
    out_file << code;
    out_file.flush();
    out_file.close();
    auto codewrite_t = t.Passed();
    LOG(INFO) << "Codegen Time: " << codewrite_t + codegen_t << "s";

    // compile and run
    auto compile_cmd = fmt::format("cmake --build {compile_path} --target prof_runner",
                                   fmt::arg("compile_path", PROJECT_BINARY_DIR));
    LOG(INFO) << "CMD: " << compile_cmd;
    t.Reset();
    int flag = system(compile_cmd.c_str());
    int num_try = 3;
    while (num_try > 0 && flag != 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        flag = system(compile_cmd.c_str());
        num_try--;
    }
    if (flag != 0) exit(-1 && "compilation error");
    auto compile_t = t.Passed();
    LOG(INFO) << "Compilation Time: " << compile_t << "s";
    auto mkdir_cmd = fmt::format("mkdir -p {bin_dir}",
                                 fmt::arg("bin_dir", PROJECT_PLAN_DIR));
    LOG(INFO) << "CMD: " << mkdir_cmd;
    flag = system(mkdir_cmd.c_str());
    if (flag != 0) exit(-1 && "mkdir error");

    std::filesystem::path bin_path = std::filesystem::path(CMAKE_RUNTIME_OUTPUT_DIRECTORY) / "prof_runner";
    std::filesystem::path dst_path = std::filesystem::path(PROJECT_PLAN_DIR) / std::to_string(config.exp_id);
    auto mv_cmd = fmt::format("mv {bin_path} {dst_path}",
                              fmt::arg("bin_path", bin_path.string()),
                              fmt::arg("dst_path", dst_path.string()));

    flag = system(mv_cmd.c_str());
    if (flag != 0) exit(-1 && "mv error");

    log.expId = config.exp_id;
    log.patternSize = sqrt(config.pat.size());
    log.compileTime = compile_t;
    log.codegenTime = codegen_t;
    log.pruningType = config.codegen.pruningType;
    log.parallelType = config.codegen.parType;
    log.adjMatType = config.codegen.adjMatType;
    log.patternAdj = config.pat;
    log.patternName = config.pattern_name;
    log.dataName = config.data_name;
    log.save(PROJECT_LOG_DIR);
};

std::string get_prefix(AppConfig config) {
    std::string adjtype = "";
    switch (config.codegen.adjMatType) {
        case minigraph::AdjMatType::VertexInduced:
            adjtype = "vertex";
            break;

        case minigraph::AdjMatType::EdgeInduced:
            adjtype = "edge";
            break;

        case minigraph::AdjMatType::EdgeInducedIEP:
            adjtype = "edgeiep";
            break;
    }
    std::string prunetype = "";
    switch (config.codegen.pruningType) {
        case minigraph::PruningType::None:
            prunetype = "none";
            break;

        case minigraph::PruningType::Static:
            prunetype = "lazy";
            break;

        case minigraph::PruningType::Eager:
            prunetype = "eager";
            break;

        case minigraph::PruningType::Online:
            prunetype = "online";
            break;

            case minigraph::PruningType::CostModel:
            prunetype = "costmodel";
            break;
    }
    std::string prefix = fmt::format("{pattern}_{graph}_{adjtype}_{prunetype}",
                                       fmt::arg("pattern", config.pattern_name),
                                       fmt::arg("graph", config.data_name),
                                       fmt::arg("adjtype", adjtype),
                                       fmt::arg("prunetype", prunetype));
    return prefix;
}

std::string get_aggr_path(AppConfig config) {
    std::filesystem::path profDir = PROJECT_PROFILE_DIR;
    std::string prefix = get_prefix(config);
    std::string filename = fmt::format("{}_aggr", prefix);
    return profDir / filename;
};

std::string get_loop_path(AppConfig config) {
    std::filesystem::path profDir = PROJECT_PROFILE_DIR;
    std::string prefix = get_prefix(config);
    std::string filename = fmt::format("{}_loop", prefix);
    return profDir / filename;
}

void run(AppConfig config) {
    std::filesystem::path bin_path = std::filesystem::path(PROJECT_PLAN_DIR) / std::to_string(config.exp_id);
    auto run_cmd = fmt::format("{bin_path} {exp_id} {data_dir} {aggr_path} {loop_path}",
                               fmt::arg("bin_path", bin_path.string()),
                               fmt::arg("data_dir", get_data_dir(config).string()),
                               fmt::arg("exp_id", config.exp_id),
                               fmt::arg("aggr_path", get_aggr_path(config)),
                               fmt::arg("loop_path", get_loop_path(config)));
    LOG(INFO) << "CMD: " << run_cmd;
    std::string run_results = exec(run_cmd.c_str());
    LOG(MSG) << run_results;
}

void compile_and_run(AppConfig config) {
    LOG(MSG) << "Experiment Setting: " << get_prefix(config);
    compile(config);
    run(config);
}

std::vector<CodeGenConfig> gen_configs() {
    std::vector<CodeGenConfig> out;
//    PruningType AllPrunningType[] = {PruningType::None, PruningType::Static, PruningType::Eager, PruningType::Online, PruningType::CostModel};
    PruningType AllPrunningType[] = {PruningType::CostModel, PruningType::Eager, PruningType::Static, PruningType::None};
    ParallelType AllParallelType[] = {ParallelType::NestedRt};
//    AdjMatType AllAdjType[] = {AdjMatType::VertexInduced, AdjMatType::EdgeInduced, AdjMatType::EdgeInducedIEP};
     AdjMatType AllAdjType[] = {AdjMatType::VertexInduced};
    for (auto adjMatType: AllAdjType) {
        for (auto parallelType: AllParallelType) {
            for (auto prunningType: AllPrunningType) {
                CodeGenConfig conf;
                conf.pruningType = prunningType;
                conf.adjMatType = adjMatType;
                conf.parType = parallelType;
                conf.runnerType = minigraph::RunnerType::Profiling;
                out.push_back(conf);
            }
        }
    }
    return out;
}

int main(int argc, char *argv[]) {
    std::string P1 = "0111101111011110";
    std::string P2 = "0110010111110110110001100";
    std::string P3 = "0110010111110110110101110";
    std::string P4 = "0111110111110111110111110";
    std::string P5 = "011111101111110110111000111000110000";
    std::string P6 = "011110101101110011110000101000011000";
    std::string P7 = "011110101011110010100001111000010100";
    std::string P8 = "0111111101111111011111110111111101111111001111100";
    std::string triangle = "011101110";
    std::string fourClique = "0111101111011110";
    std::string fiveClique = "0111110111110111110111110";

    std::vector<std::pair<std::string, std::string>> patterns = {
            {"p1", P1},{"p2", P2},
            {"p3", P3},{"p4", P4},
            {"p5", P5},{"p6", P6},
            {"p7", P7},{"p8", P8}
    };
    std::vector<std::string> small_graphs = {"wiki", "patents", "youtube"};
    std::vector<std::string> large_graphs = {"lj", "orkut", "friendster"};
    std::vector<std::pair<std::string, std::string>> small_patterns = {{"p1", P1},
                                                                       {"p2", P2},
                                                                       {"p3", P3},
                                                                       {"p4", P4}};
    std::vector<std::pair<std::string, std::string>> large_patterns = {{"p5", P5},
                                                                       {"p6", P6},
                                                                       {"p7", P7},
                                                                       {"p8", P8}};

    std::vector<std::pair<std::string, std::string>> foo_patterns = {{"p3", P3},
                                                                    {"p4", P4},
                                                                    {"p5", P5}};
    std::vector<AppConfig> configs;
    auto code_gen_configs = gen_configs();
    int exp_id = 0;

    std::vector<std::string> debug_graphs = {"lj"};
    // std::vector<std::pair<std::string, std::string>> debug_patterns = {{"3clique", triangle}, {"4clique", fourClique}, {"5clique", fiveClique}};
    std::vector<std::pair<std::string, std::string>> debug_patterns = {{"p7", P7}};

    for (auto graph: debug_graphs) {
        for (auto pattern: debug_patterns) {
            for (auto code_gen: code_gen_configs) {
                AppConfig config;
                config.exp_id = exp_id++;
                config.pattern_name = pattern.first;
                config.pat = pattern.second;
                config.codegen = code_gen;
                config.data_name = graph;
                configs.push_back(config);
            }
        }
    }
//
//    for (auto graph: large_graphs) {
//        for (auto pattern: patterns) {
//            for (auto code_gen: code_gen_configs) {
//                AppConfig config;
//                config.exp_id = exp_id++;
//                config.pattern_name = pattern.first;
//                config.pat = pattern.second;
//                config.codegen = code_gen;
//                config.data_name = graph;
//                configs.push_back(config);
//            }
//        }
//    }
    for (auto config: configs) {
        compile_and_run(config);
    }
}
