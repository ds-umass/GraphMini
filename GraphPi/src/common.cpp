#include "common.h"
#include <sys/time.h>
#include <cstdlib>

double get_wall_time() {
    struct timeval time;
    if(gettimeofday(&time,NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

void PatternType_printer(PatternType type) {
    if(type == PatternType::Rectangle) {
        printf("Rectangle\n");
    }
    if(type == PatternType::Pentagon) {
        printf("Pentagon\n");
    }
    if(type == PatternType::House) {
        printf("House\n");
    }
    if(type == PatternType::Hourglass) {
        printf("Hourglass\n");
    }
    if(type == PatternType::Cycle_6_Tri) {
        printf("Cycle_6_Tri\n");
    }
    if(type == PatternType::Clique_7_Minus) {
        printf("Clique_7_Minus\n");
    }
}

bool is_equal_adj_mat(const int* adj_mat1, const int* adj_mat2, int size) {
    for(int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            if(adj_mat1[INDEX(i,j,size)] != adj_mat2[INDEX(i,j,size)]) 
                return false;
    return true;
}

// modified method
void GetDataType(DataType &type, const std::string str) {
    type = DataType::Invalid;
    if(str == "Wiki-Vote" || str == "wiki") {
        type = DataType::Wiki_Vote;
    }

    if( str == "Patents" || str == "patents") {
        type = DataType::Patents;
    }

    if (str == "YouTube" || str == "youtube") {
        type = DataType::YouTube;
    }

    if(str == "Orkut" || str == "orkut") {
        type = DataType::Orkut;
    }

    if(str == "LiveJournal" || str == "lj") {
        type = DataType::LiveJournal;
    }

    if (str == "Friendster" || str == "friendster") {
        type = DataType::Friendster;
    }

    if(str == "MiCo" || str == "mico") {
        type = DataType::MiCo;
    }
    if(str == "Twitter" || str == "twitter") {
        type = DataType::Twitter;
    }
    if(str == "CiteSeer" || str == "citeseer" ) {
        type = DataType::CiteSeer;
    }

    if(str == "complete8" ) {
        type = DataType::complete8;
    }
}

int read_int() {
    char ch = getchar();
    while((ch < '0' || ch > '9') && ch !='-') ch = getchar();
    int tag = 1;
    if(ch == '-') tag = -1, ch = getchar();
    int x = 0;
    while( ch >= '0' && ch <= '9') x = x* 10 + ch -'0', ch = getchar();
    return x * tag;
}

unsigned int read_unsigned_int() {
    char ch = getchar();
    while((ch < '0' || ch > '9') ) ch = getchar();
    unsigned int x = 0;
    while( ch >= '0' && ch <= '9') x = x* 10 + ch -'0', ch = getchar();
    return x;
}
