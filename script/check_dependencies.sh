#!/bin/bash
# 依赖检查脚本

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}=== 依赖检查 ===${NC}"
echo ""

# 检查函数
check_command() {
    local cmd=$1
    local name=$2
    local install_cmd=$3
    
    if command -v $cmd &> /dev/null; then
        local version=$($cmd --version 2>&1 | head -n 1)
        echo -e "${GREEN}✓${NC} $name: $version"
        return 0
    else
        echo -e "${RED}✗${NC} $name: 未安装"
        if [ ! -z "$install_cmd" ]; then
            echo -e "  ${YELLOW}安装命令:${NC} $install_cmd"
        fi
        return 1
    fi
}

check_pkg_config() {
    local pkg=$1
    local name=$2
    local install_cmd=$3
    
    if pkg-config --exists $pkg 2>/dev/null; then
        local version=$(pkg-config --modversion $pkg 2>/dev/null)
        echo -e "${GREEN}✓${NC} $name: $version"
        return 0
    else
        echo -e "${RED}✗${NC} $name: 未安装"
        if [ ! -z "$install_cmd" ]; then
            echo -e "  ${YELLOW}安装命令:${NC} $install_cmd"
        fi
        return 1
    fi
}

# 检查编译工具
echo -e "${BLUE}编译工具:${NC}"
check_command cmake "CMake" "sudo apt install cmake"
check_command g++ "GCC" "sudo apt install g++"
check_command make "Make" "sudo apt install make"
echo ""

# 检查库
echo -e "${BLUE}必需库:${NC}"
check_pkg_config eigen3 "Eigen3" "sudo apt install libeigen3-dev"
check_pkg_config ceres "Ceres Solver" "sudo apt install libceres-dev"
echo ""

# 检查可选工具
echo -e "${BLUE}可选工具:${NC}"
check_command git "Git" "sudo apt install git"
check_command gdb "GDB" "sudo apt install gdb"
echo ""

# 总结
echo -e "${BLUE}=== 检查完成 ===${NC}"
echo ""
echo -e "${YELLOW}提示:${NC}"
echo "  - 如果有红色的✗，请运行相应的安装命令"
echo "  - 在WSL2 Ubuntu终端中运行安装命令"
echo "  - 安装后可能需要重新加载终端环境"
echo ""
