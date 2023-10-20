cmake_minimum_required(VERSION 3.24)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(Helm)

include(FetchContent)

FetchContent_Declare(
  eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG e7248b26a1ed53fa030c5c459f7ea095dfd276ac
  )

FetchContent_MakeAvailable(eigen3)

find_package(OpenMP REQUIRED)

aux_source_directory("${PROJECT_SOURCE_DIR}/src" Helm_SOURCE)
add_executable(Helm ${Helm_SOURCE})
target_link_libraries(Helm PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)


