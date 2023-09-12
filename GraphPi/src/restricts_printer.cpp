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
#include <vector>

int main(int argc,char *argv[]) {
   
    int size = atoi(argv[1]);
    char* adj_mat = argv[2];

    Pattern pattern(size, adj_mat);
    bool is_valid;
    Schedule schedule(pattern, is_valid, 0, 0, 0, 0, 0, 0);
    std::vector< std::vector< std::pair<int,int> > > restricts;
    schedule.restricts_generate(schedule.get_adj_mat_ptr(), restricts);
    schedule.print_schedule();
    for(size_t i = 0; i < restricts.size(); i++) {
        auto& r = restricts.at(i); 
        for(size_t j = 0; j < r.size(); j++) {
            auto& p = r.at(j);
            printf("(%d,%d)",p.first,p.second);
        }
        puts("");
    }
    return 0;
}
