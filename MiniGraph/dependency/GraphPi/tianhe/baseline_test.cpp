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
#include <math.h>

void test_pattern(Graph *g, const Pattern &pattern, int performance_modeling_type, int restricts_type,
                  bool use_in_exclusion_optimize = false) {
    int thread_num = 32;
    double t1, t2;

    bool is_pattern_valid;
    Schedule schedule(pattern, is_pattern_valid, performance_modeling_type, restricts_type, use_in_exclusion_optimize,
                      g->v_cnt, g->e_cnt, g->tri_cnt);
    assert(is_pattern_valid);

    t1 = get_wall_time();
    long long ans = g->pattern_matching(schedule, thread_num);
    t2 = get_wall_time();

    printf("ans %lld\n", ans);
    printf("time %.6lf\n", t2 - t1);
    schedule.print_schedule();
    const auto &pairs = schedule.restrict_pair;
    printf("%d ", pairs.size());
    for (auto &p: pairs)
        printf("(%d,%d)", p.first, p.second);
    puts("");
    fflush(stdout);
}

std::string p1 = "0111101011011010";
std::string p2 = "011110101101110011110000101000011000";
std::string p3 = "011111101111110110111000111000110000";
std::string p4 = "011110101011110010100001111000010100";
std::string p5 = "0111111101111111011001110110111100011010001100000";

int main(int argc, char *argv[]) {

//    if (argc < 5) {
//        printf("Usage: %s dataset_name graph_file pattern_size pattern_adjacency_matrix\n", argv[0]);
//        printf("Example(Triangle counting on dataset WikiVote) : \n");
//        printf("%s Wiki-Vote ../../dataset/wiki-vote_input 3 011101110\n", argv[0]);
//        return 0;
//    }

    Graph *g;
    DataLoader D;

    const std::string type = "Wiki-Vote";
    const std::string path = "/home/ubuntu/MiniGraph/dependency/GraphPi/dataset/wiki-vote_input";

    int size = int(sqrt(p2.size()));
    char *adj_mat = const_cast<char *>(p2.c_str());

    // comments in include/schedule.h explain the meaning of these parameters.
    int test_type = 1; // performance_modeling_type = restricts_type = use_in_exclusion_optimize = 1

    DataType my_type;

    GetDataType(my_type, type);

    if (my_type == DataType::Invalid) {
        printf("Dataset not found!\n");
        return 0;
    }

    assert(D.load_data(g, my_type, path.c_str()) == true);

    printf("Load data success!\n");
    fflush(stdout);

    Pattern p(size, adj_mat);
    test_pattern(g, p, test_type, test_type, 0);

    delete g;
}
