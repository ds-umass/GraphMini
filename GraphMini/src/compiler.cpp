//
// Created by ubuntu on 1/1/23.
//
#include "codegen.h"
#include "configure.h"
#include "common.h"
#include <vector>
#include <iostream>
#include <fmt/format.h>
#include <filesystem>
#include <fstream>
using namespace minigraph;


//std::string adjmat_to_graphzero_pattern_list(std::string pat){
//    std::ostringstream header, adj_list;
//    int p_size = (int) sqrt(pat.size());
//    if(p_size * p_size != pat.size()) exit(-1 && "Invalid size");
//    int num_edges = 0;
//    for (size_t i = 0; i < p_size; i++){
//        for (size_t j = 0; j < p_size; j++){
//            if (i == j) continue;
//            if (pat.at(i * p_size + j) == '1') {
//                num_edges++;
//                adj_list << i << " " << j << "\n";
//            }
//        }
//    }
//
//    header << "1\n" << p_size << " " << num_edges << "\n";
//    for (int i = 0; i < p_size - 1; i++) {
//        header << "0 ";
//    }
//    header << "0\n";
//    return header.str() + adj_list.str();
//}
struct AppConfig
{
    CodeGenConfig code_gen_config;
    std::filesystem::path data_dir, output_dir;
    std::string pattern_name;
};

void run_mg(std::string adj_mat,
            std::filesystem::path data_dir,
            CodeGenConfig config){
    MetaData meta;
    meta.read(data_dir);
    std::string code = gen_code(adj_mat, config, meta);
    std::filesystem::path code_file(PROJECT_SOURCE_DIR);
    code_file /= "src";
    code_file /= "codegen_output";
    code_file /= "generated_plan.cpp";
    std::ofstream out_file(code_file);
    out_file << code;
    out_file.close();

//    std::cout << "----------GENERATED CODE BEGIN---------\n";
//    std::cout << code;
//    std::cout << "----------GENERATED CODE END---------\n";

    // compile and run
    Timer t;
    std::filesystem::path compile_dir = PROJECT_BINARY_DIR;
    int flag = system(fmt::format("cmake --build {compile_path} --target runner",
                                  fmt::arg("compile_path", compile_dir.string())).c_str());
    if (flag != 0) exit(-1 && "compilation error");
    LOG(INFO) << "Binary Compilation Time: " << t.Passed() << "s";

    std::filesystem::path bin_path = std::filesystem::path(CMAKE_RUNTIME_OUTPUT_DIRECTORY) / "runner";
    flag = system(fmt::format("{bin_path} 0 {data_dir}",
                       fmt::arg("bin_path", bin_path.string()),
                       fmt::arg("data_dir", data_dir.string())).c_str());
    if (flag != 0) exit(flag && "codegen_output error");
};

int main(int argc, char * argv[]){
    std::string p0 = "0111101011011010";
    std::string p1 = "0111110111110111110111110";
    std::string p2 = "011110101101110011110000101000011000";
    std::string p3 = "011111101111110110111000111000110000";
    std::string p4 = "011110101011110010100001111000010100";
    std::string p5 = "0111111101111111011001110110111100011010001100000";
    std::string p6 = "0111100101110011011001110111111101100011000001100";
    std::string P5 = "0111110111110111110111110";
    std::string P6 = "011010101001110100001011100101010110";

    std::filesystem::path data_dir("/home/ubuntu");
    data_dir /= "dataset";
    data_dir /= "wiki";
    CodeGenConfig config;
    config.adjMatType = minigraph::AdjMatType::VertexInduced;
    config.parType = minigraph::ParallelType::OpenMP;
    config.pruningType = minigraph::PruningType::None;
    run_mg(P6, data_dir, config);
}