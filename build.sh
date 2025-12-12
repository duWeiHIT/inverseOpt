#!/bin/bash
# 构建脚本 - 用于WSL2 Ubuntu环境

set -e

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}=== 热传导反演求解器 - 构建脚本 ===${NC}"
echo ""

# 检查CMake
if ! command -v cmake &> /dev/null; then
    echo -e "${RED}错误: 未找到CMake${NC}"
    echo "请运行: sudo apt install cmake"
    exit 1
fi

# 检查编译器
if ! command -v g++ &> /dev/null; then
    echo -e "${RED}错误: 未找到g++编译器${NC}"
    echo "请运行: sudo apt install g++"
    exit 1
fi

# 检查Eigen3
if ! pkg-config --exists eigen3 2>/dev/null && [ ! -d "/usr/include/eigen3" ]; then
    echo -e "${YELLOW}警告: 未找到Eigen3${NC}"
    echo "请运行: sudo apt install libeigen3-dev"
else
    echo -e "${GREEN}✓${NC} Eigen3: 已安装"
fi

# 检查Ceres (使用多种方法检查)
CERES_FOUND=false
if pkg-config --exists ceres 2>/dev/null; then
    CERES_FOUND=true
    echo -e "${GREEN}✓${NC} Ceres Solver: $(pkg-config --modversion ceres 2>/dev/null)"
elif [ -f "/usr/lib/cmake/Ceres/CeresConfig.cmake" ] || [ -f "/usr/lib/x86_64-linux-gnu/cmake/Ceres/CeresConfig.cmake" ]; then
    CERES_FOUND=true
    echo -e "${GREEN}✓${NC} Ceres Solver: 已安装 (通过CMake配置检测)"
elif ldconfig -p | grep -q libceres 2>/dev/null; then
    CERES_FOUND=true
    echo -e "${GREEN}✓${NC} Ceres Solver: 已安装 (通过库文件检测)"
fi

if [ "$CERES_FOUND" = false ]; then
    echo -e "${YELLOW}警告: 未找到Ceres Solver${NC}"
    echo "请运行: sudo apt install libceres-dev"
fi

echo ""
echo -e "${YELLOW}CMake版本:${NC}"
cmake --version | head -n 1

echo ""
echo -e "${YELLOW}编译器版本:${NC}"
g++ --version | head -n 1

echo ""
echo -e "${YELLOW}正在创建构建目录...${NC}"
mkdir -p build
cd build

echo -e "${YELLOW}正在配置项目...${NC}"
if cmake -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-O3 -march=native" \
    .. ; then
    echo -e "${GREEN}✓${NC} CMake配置成功"
else
    echo -e "${RED}CMake配置失败${NC}"
    exit 1
fi

echo ""
echo -e "${YELLOW}正在编译项目...${NC}"
if cmake --build . --config Release; then
    echo -e "${GREEN}✓${NC} 编译成功"
else
    echo -e "${RED}编译失败${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}=== 编译完成 ===${NC}"
echo -e "${GREEN}可执行程序位置:${NC}"
echo "  - $(pwd)/bin/singleLayerInverse"
echo "  - $(pwd)/bin/multiLayerInverse"
echo "  - $(pwd)/bin/dataGen"
echo ""
