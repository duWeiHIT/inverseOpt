cmake_minimum_required(VERSION 3.10)
project(ThermalInversion CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 查找必需的包
find_package(Eigen3 REQUIRED NO_MODULE)

# 优先使用apt安装的Ceres
set(CMAKE_PREFIX_PATH "/usr/lib/cmake/Ceres")
set(Ceres_DIR "/usr/lib/cmake/Ceres")
find_package(Ceres REQUIRED)

# 如果找不到，尝试pkg-config
if(NOT Ceres_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CERES REQUIRED ceres)
endif()

# 包含目录
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# 为每个源文件创建可执行程序
add_executable(singleLayerInverse src/singleLayerInverse.cpp)
target_link_libraries(singleLayerInverse PRIVATE Eigen3::Eigen ceres)

add_executable(multiLayerInverse src/multiLayerInverse.cpp)
target_link_libraries(multiLayerInverse PRIVATE Eigen3::Eigen ceres)


# 编译选项
if(MSVC)
    target_compile_options(singleLayerInverse PRIVATE /W4)
    target_compile_options(multiLayerInverse PRIVATE /W4)
else()
    target_compile_options(singleLayerInverse PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(multiLayerInverse PRIVATE -Wall -Wextra -Wpedantic)
endif()
