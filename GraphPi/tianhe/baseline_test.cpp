#include <../include/graph.h>
#include <../include/dataloader.h>
#include "../include/pattern.h"
#include "../include/schedule.h"
#include "../include/common.h"
#include "../include/motif_generator.h"
#include <omp.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <thread>

void test_pattern(Graph* g, const Pattern &pattern, int performance_modeling_type, int restricts_type, bool use_in_exclusion_optimize = false) {
    int thread_num = std::thread::hardware_concurrency();
    double t0,t1,t2;
    bool is_pattern_valid;

    t0 = get_wall_time();
    Schedule schedule(pattern, is_pattern_valid, performance_modeling_type, restricts_type, use_in_exclusion_optimize, g->v_cnt, g->e_cnt, g->tri_cnt);
    assert(is_pattern_valid);

    t1 = get_wall_time();
    long long ans = g->pattern_matching(schedule, thread_num);
    t2 = get_wall_time();

    schedule.print_schedule();
    printf("Restricts: ");
    const auto& pairs = schedule.restrict_pair;
    printf("%ld ",pairs.size());
    for(auto& p : pairs)
        printf("(%d,%d)",p.first,p.second);
    puts("");
    fflush(stdout);

    printf("RESULT=%lld\n", ans);
    printf("CODE_GENERATION_TIME(s)=%.6lf\n", t1 - t0);
    printf("CODE_EXECUTION_TIME(s)=%.6lf\n", t2 - t1);
    printf("NUM_THREADS=%d\n",thread_num);
}

int main(int argc,char *argv[]) {
    if(argc < 6) {
        printf("Usage: %s dataset_name graph_file pattern_size pattern_adjacency_matrix\n", argv[0]);
        printf("Example(Triangle counting on dataset WikiVote) : \n");
        printf("%s Wiki-Vote ../../dataset/wiki-vote_input 3 011101110 [1 if USE_IEP else 0]\n", argv[0]);
        return 0;
    }

    Graph *g;
    DataLoader D;

    const std::string type = argv[1];
    const std::string path = argv[2];
    
    int size = atoi(argv[3]);
    char* adj_mat = argv[4];

    // comments in include/schedule.h explain the meaning of these parameters.
    int performance_modeling_type = 1;
    int restricts_type = 1;
    int use_in_exclusion_optimize = std::atoi(argv[5]);
    DataType my_type;
    
    GetDataType(my_type, type);

    if(my_type == DataType::Invalid) {
        printf("Dataset not found!\n");
        return 0;
    }

    assert(D.load_data(g,my_type,path.c_str())==true); 

    printf("Load data success!\n");
    fflush(stdout);

    Pattern p(size, adj_mat);
    test_pattern(g, p, performance_modeling_type, restricts_type, use_in_exclusion_optimize);
    
    delete g;
}
