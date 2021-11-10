set(ENGINE_FILE_LIST
        matchingcommand.cpp matchingcommand.h
        engine.cpp
        execution_tree.h
        execution_tree_node.h
        execution_tree_node.cpp
        execution_tree_generator.cpp
        execution_tree_generator.h
        execution_tree.cpp
        preprocessor.cpp
        preprocessor.h
        encoder.cpp
        encoder.h
        query_plan_generator.cpp
        query_plan_generator.h)

add_executable(RapidMatch.out ${ENGINE_FILE_LIST})
add_executable(Playground.out playground.cpp)
target_link_libraries(RapidMatch.out graph utility nd)