project(testCircleDecomposition)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})

add_definitions(-g)

LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
TARGET_LINK_LIBRARIES(testCircleDecomposition ${DGTAL_LIBRARIES})
