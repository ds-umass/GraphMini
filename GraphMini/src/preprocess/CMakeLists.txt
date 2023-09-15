add_library(preprocess STATIC
        graph_converter.cpp)

target_link_libraries(preprocess PRIVATE common TBB::tbb TBB::tbbmalloc)