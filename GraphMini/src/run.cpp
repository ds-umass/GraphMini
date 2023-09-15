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
using namespace minigraph;

std::string exec(const char *cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}


struct AppConfig {
    int exp_id;
    CodeGenConfig codegen;
    std::string data_name;
    std::string pattern_name;
    std::string pat;
    std::string graph_dir;
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
    } else if (config.codegen.pruningType == PruningType::CostModel) {
        out /= prefix + "costmodel.txt";
    }
    return out;
}

std::filesystem::path code_path() {
    std::filesystem::path code_file(PROJECT_SOURCE_DIR);
    code_file /= "src";
    code_file /= "codegen_output";
    code_file /= "plan.cpp";
    return code_file;
}

void compile(AppConfig config) {
    CompilerLog log;
    MetaData meta;
    meta.read(config.graph_dir);
    Timer t;
    std::string code = gen_code(config.pat, config.codegen, meta);
    auto codegen_t = t.Passed();
    t.Reset();
    std::ofstream out_file(code_path());

    LOG(MSG) << "Code is written to: " << code_path();
    out_file << code;
    out_file.flush();
    out_file.close();
    auto codewrite_t = t.Passed();

    // format code
    auto reformat_cmd = fmt::format("clang-format -i {}", code_path().string());
    auto reformat_flag = system(reformat_cmd.c_str());
    if (reformat_flag != 0) {
        LOG(MSG) << "Install clang-format to have formatted code at: " << code_path();
    }

    LOG(MSG) << "CODE_GENERATION_TIME(s)=" << codewrite_t + codegen_t;

    // compile and run
    auto compile_cmd = fmt::format("cmake --build {compile_path} --target runner 1>>/dev/null 2>>/dev/null",
                                   fmt::arg("compile_path", PROJECT_BINARY_DIR));
    // LOG(MSG) << "CMD: " << compile_cmd;
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
    LOG(MSG) << "COMPILATION_TIME(s)=" << compile_t;
    auto mkdir_cmd = fmt::format("mkdir -p {bin_dir}",
                                 fmt::arg("bin_dir", PROJECT_PLAN_DIR));
    // LOG(MSG) << "CMD: " << mkdir_cmd;
    flag = system(mkdir_cmd.c_str());
    if (flag != 0) exit(-1 && "mkdir error");

    std::filesystem::path bin_path = std::filesystem::path(CMAKE_RUNTIME_OUTPUT_DIRECTORY) / "runner";
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

void run(AppConfig config) {
    std::filesystem::path bin_path = std::filesystem::path(PROJECT_PLAN_DIR) / std::to_string(config.exp_id);
    auto run_cmd = fmt::format("{bin_path} {exp_id} {data_dir}",
                               fmt::arg("bin_path", bin_path.string()),
                               fmt::arg("data_dir", config.graph_dir),
                               fmt::arg("exp_id", config.exp_id));
    std::string run_results = exec(run_cmd.c_str());
    LOG(MSG) << run_results;
}

void compile_and_run(AppConfig config) {
    compile(config);
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(500ms);
    run(config);
}

int main(int argc, char *argv[]) {
    using namespace minigraph;
    if (argc < 8) {
        std::cout << "./run [graph_name] [graph_dir] [query_name] [query] [adj_type] [prun_type] [par_type] [exp_id=-1 (optional)]\n";
        std::cout << "adj_type: 0=VertexInduced; 1=EdgeInduced; 2=EdgeInducedIEP\n";
        std::cout << "prun_type: 0=None; 1=Static; 2=Eager; 3=Online; 4=CostModel\n";
        std::cout << "par_type: 0=OpenMP; 1=TbbTop; 2=Nested; 3=NestedRt\n";
        std::cout << "For example:\n./MiniGraph/build/bin/run wiki ./Datasets/MiniGraph/wiki/ P1 0111101111011110 0 4 3\n";
        return 0;
    }

    std::string graph_name={argv[1]};
    std::string graph_dir={argv[2]};
    std::string query_name={argv[3]};
    std::string query_str={argv[4]};
    int adjmat_type_int = std::atoi(argv[5]);
    int prun_type_int = std::atoi(argv[6]);
    int par_type_int = std::atoi(argv[7]);
    int exp_id = -1;
    if (argc == 9) exp_id = std::atoi(argv[8]);

    AdjMatType adjmat_type = (AdjMatType) (adjmat_type_int);
    PruningType prun_type = static_cast<PruningType>(prun_type_int);
    ParallelType par_type = static_cast<ParallelType>(par_type_int);

    CodeGenConfig conf;
    conf.adjMatType  = adjmat_type;
    conf.pruningType = prun_type;
    conf.parType     = par_type;

    AppConfig config;
    config.exp_id = exp_id;
    config.pattern_name = query_name;
    config.pat = query_str;
    config.codegen = conf;
    config.data_name = graph_name;
    config.graph_dir = graph_dir;
    compile_and_run(config);
}
