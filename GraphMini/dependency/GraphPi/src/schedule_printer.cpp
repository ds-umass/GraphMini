#include <../include/graph.h>
#include <../include/dataloader.h>
#include "../include/pattern.h"
#include "../include/schedule.h"
#include "../include/common.h"
#include "../include/motif_generator.h"

#include <assert.h>
#include <iostream>
#include <string>
#include <algorithm>

#include <chrono>

using namespace std::chrono;
//void print_schedule(const Pattern &pattern, int v_num, int e_num, int tri_num, int performance_modeling_type = 1){
//    bool is_pattern_valid;
//    int restricts_type = performance_modeling_type;
//    int use_in_exclusion_optimize = 0;
//    Schedule schedule(pattern, is_pattern_valid, performance_modeling_type, restricts_type, use_in_exclusion_optimize, v_num, e_num, tri_num);
//    schedule.print_schedule();
//};

void print_graphpi_schedule(const char *adj_mat, int size, int v_num, int e_num, int tri_num) {
    Schedule schedule = Schedule();
    schedule.get_schedule(adj_mat, size, v_num, e_num, tri_num);
}

void print_greedy_schedule(const char *adj_mat, int size) {
    Schedule schedule = Schedule();
    schedule.get_greedy_schedule(adj_mat, size);
}

int main(int argc, char *argv[]) {

    if (argc != 7) {
        printf("./scheduler_printer [pattern_size] [adj_mat] [v_num] [e_num] [tri_num] [model (GraphPi / minigraph)]");
        return 0;
    }
    char *adj_mat;
    int pattern_size, v_num, e_num, tri_num;
    pattern_size = atoi(argv[1]);
    adj_mat = argv[2];
    v_num = atoi(argv[3]);
    e_num = atoi(argv[4]);
    tri_num = atoi(argv[5]);
    std::string model(argv[6]);
    if (model == "GraphPi") {
        print_graphpi_schedule(adj_mat, pattern_size, v_num, e_num, tri_num);
    } else if (model == "minigraph") {
        print_greedy_schedule(adj_mat, pattern_size);
    }
    return 0;
}
