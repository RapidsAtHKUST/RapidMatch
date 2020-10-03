add_executable(ReassignLabel ReassignLabel.cpp)
add_executable(SelectEdge SelectEdge.cpp)
add_executable(GraphConverter.out GraphConverter.cpp)
target_link_libraries(GraphConverter.out graph utility)

add_executable(EdgeListConverter.out EdgeListToCSR.cpp)