add_library(codegen STATIC
        codegen.cpp
        ir.cpp)

target_link_libraries(codegen PRIVATE common graph_mining fmt::fmt)
