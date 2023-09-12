//INCLUDE_HEADER_HERE
#include "support.h"
#include <vector>
#include <iostream>
#include <chrono>
#include <omp.h>

int main(int argc, char **argv) {
  /**if(argc != 2) {
    std::cerr << "usage: ./test <graph file>\n";
    exit(1);
    }*/
  // std::cout << argv[1] << "\n";
	auto t1 = std::chrono::high_resolution_clock::now();
  Graph g(argv[1]);
	auto t2 = std::chrono::high_resolution_clock::now();
	std::cout << "GRAPH_LOADING_TIME(S)=" << std::chrono::duration<double, std::ratio<1,1>>(t2-t1).count() << "\n";
  std::cout.flush();
	auto t3 = std::chrono::high_resolution_clock::now();

  // INCLUDE_HERE

	auto t4 = std::chrono::high_resolution_clock::now();
	std::cout << "CODE_EXECUTION_TIME(s)=" << std::chrono::duration<double, std::ratio<1,1>>(t4-t3).count() << "\n";
  std::cout << std::flush;
  uint64_t sum = 0;
  for(uint32_t result_idx = 0; result_idx < n_counters; result_idx++) {
    for(uint32_t thread_idx = 0; thread_idx < global_counters.size(); thread_idx++) {
      if(global_counters.at(thread_idx).size() > result_idx) {
        sum += global_counters.at(thread_idx).at(result_idx);
      }
    }
  }
  std::cout << "RESULT=" << sum << "\n";
  // std::cout << "data_complexity = " << data_complexity << "\n";
  // std::cout << "time_complexity = " << time_complexity << "\n";
}
