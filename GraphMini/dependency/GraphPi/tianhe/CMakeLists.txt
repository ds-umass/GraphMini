ADD_EXECUTABLE(baseline_test baseline_test.cpp)
TARGET_LINK_LIBRARIES(baseline_test graph_mining)

#ADD_EXECUTABLE(in_exclusion_optimize_test in_exclusion_optimize_test.cpp)
#TARGET_LINK_LIBRARIES(in_exclusion_optimize_test graph_mining)

#ADD_EXECUTABLE(all_schedule_test all_schedule_test.cpp)
#TARGET_LINK_LIBRARIES(all_schedule_test graph_mining)
#
#ADD_EXECUTABLE(performance_test performance_test.cpp)
#TARGET_LINK_LIBRARIES(performance_test graph_mining)
#
#ADD_EXECUTABLE(schedules_test schedules_test.cpp)
#TARGET_LINK_LIBRARIES(schedules_test graph_mining)

ADD_EXECUTABLE(code_generator code_generator.cpp)
TARGET_LINK_LIBRARIES(code_generator graph_mining)

#ADD_EXECUTABLE(unfold_test unfold_test.cpp)
#TARGET_LINK_LIBRARIES(unfold_test graph_mining)
