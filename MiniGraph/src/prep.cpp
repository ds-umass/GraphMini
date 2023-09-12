//
// Created by ubuntu on 1/2/23.
//

#include "common.h"
#include "preprocess/graph_converter.h"
int main(int argc, char * argv[]){
    using namespace minigraph;
    CHECK(argc == 2) << "Usage: ./prep [directory contain snap.txt downloaded from SNAP]";
    minigraph::GraphConverter converter;
    std::filesystem::path in_dir{ argv[1] };
    converter.convert(in_dir);
}