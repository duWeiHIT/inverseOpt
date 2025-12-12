#!/bin/bash
# 高级CMake配置脚本 - WSL2 Ubuntu环境

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# 默认值
BUILD_TYPE="Release"
ENABLE_TESTS=OFF
ENABLE_DOCS=OFF
PARALLEL_JOBS=$(nproc)
INSTALL_PREFIX="/usr/local"

# 显示帮助信息
show_help() {
    cat << EOF
${BLUE}热传导反演求解器 - WSL2 Ubuntu CMake配置脚本${NC}

用法: ./cmake_config.sh [选项]

选项:
    -h, --help              显示此帮助信息
    -d, --debug             调试模式编译 (默认: Release)
    -t, --tests             启用测试 (默认: OFF)
    -j, --jobs N            并行编译数 (默认: $(nproc))
    -p, --prefix PATH       安装前缀 (默认: /usr/local)
    --clean                 清理旧的构建目录
    --verbose               详细输出

示例:
    ./cmake_config.sh                    # 默认Release构建
    ./cmake_config.sh -d                 # Debug构建
    ./cmake_config.sh -d -j 4            # Debug构建，4个并行任务
    ./cmake_config.sh --clean -d         # 清理后Debug构建

EOF
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -d|--debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        -t|--tests)
            ENABLE_TESTS=ON
            shift
            ;;
        -j|--jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        -p|--prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        --clean)
            echo -e "${YELLOW}清理旧的构建目录...${NC}"
            rm -rf build
            shift
            ;;
        --verbose)
            CMAKE_VERBOSE="--verbose"
            shift
            ;;
        *)
            echo -e "${RED}未知选项: $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# 显示配置信息
echo -e "${BLUE}=== 构建配置 ===${NC}"
echo -e "${GREEN}构建类型:${NC} $BUILD_TYPE"
echo -e "${GREEN}并行任务:${NC} $PARALLEL_JOBS"
echo -e "${GREEN}安装前缀:${NC} $INSTALL_PREFIX"
echo -e "${GREEN}启用测试:${NC} $ENABLE_TESTS"
echo ""

# 创建构建目录
mkdir -p build
cd build

# 运行CMake
echo -e "${YELLOW}正在配置项目...${NC}"
cmake \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
    -DENABLE_TESTS=$ENABLE_TESTS \
    -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -DNDEBUG" \
    -DCMAKE_CXX_FLAGS_DEBUG="-g -O0 -Wall -Wextra" \
    ..

if [ $? -ne 0 ]; then
    echo -e "${RED}CMake配置失败${NC}"
    exit 1
fi

# 编译
echo ""
echo -e "${YELLOW}正在编译项目 (使用 $PARALLEL_JOBS 个并行任务)...${NC}"
cmake --build . --config $BUILD_TYPE --parallel $PARALLEL_JOBS $CMAKE_VERBOSE

if [ $? -ne 0 ]; then
    echo -e "${RED}编译失败${NC}"
    exit 1
fi

# 显示完成信息
echo ""
echo -e "${GREEN}=== 编译完成 ===${NC}"
echo -e "${GREEN}可执行程序位置:${NC}"
echo "  - $(pwd)/bin/singleLayerInverse"
echo "  - $(pwd)/bin/multiLayerInverse"
echo "  - $(pwd)/bin/dataGen"
echo ""

# 显示后续步骤
echo -e "${YELLOW}后续步骤:${NC}"
echo "  1. 运行程序: ./bin/singleLayerInverse"
echo "  2. 安装: cmake --install ."
echo "  3. 清理: cd .. && rm -rf build"
echo ""
