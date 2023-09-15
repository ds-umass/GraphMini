add_library(plan STATIC
        plan.cpp plan.h)
target_precompile_headers(plan PUBLIC ../backend/backend.h)
target_link_libraries(plan PUBLIC OpenMP::OpenMP_CXX TBB::tbb TBB::tbbmalloc)

add_library(plan_profile STATIC
        plan_profile.cpp plan_profile.h)
target_precompile_headers(plan_profile PUBLIC ../backend_prof/backend.h)
target_link_libraries(plan_profile PUBLIC OpenMP::OpenMP_CXX TBB::tbb TBB::tbbmalloc)
