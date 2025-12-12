# 热传导反演求解器 - 构建指南

## 项目概述
这是一个基于C++的一维非稳态热传导反演求解器，支持多种边界条件和反演算法。

## 系统要求

### 必需软件

### 需要在Linux中安装如下依赖
- **CMake** >= 3.10
- **C++编译器** (支持C++17)
- **Eigen3** >= 3.3 (线性代数库)
- **Ceres Solver** (非线性优化库)
```bash

# 更新包管理器
sudo apt update

# 安装编译工具
sudo apt install gcc g++ cmake make

# 安装必需库
sudo apt install libeigen3-dev libceres-dev
```

## 构建步骤

### 使用Bash脚本

在WSL2 Ubuntu终端中运行:

```bash
cd /path/to/project
chmod +x build.sh
./build.sh
```

## 输出文件

编译成功后，可执行程序位于 `build/bin/` 目录:

- `singleLayerInverse` - 单层材料反演程序
- `multiLayerInverse` - 多层材料反演程序
- `dataGen` - 数据生成工具

## 运行程序

```bash
cd build/bin
./singleLayerInverse
```

## 项目结构

```
Opt/
├── CMakeLists.txt              # CMake配置文件
├── build.sh                    # Ubuntu构建脚本
├── cmake_config.sh             # 高级配置脚本
├── README.md                   # 本文件
├── src/                        # 源代码目录
│   ├── singleLayerInverse.cpp      # 单层反演主程序
│   ├── multiLayerInverse.cpp       # 多层反演主程序
│   ├── dataGen.cpp                 # 数据生成工具
│   ├── singleLayerSolver.h         # 单层求解器头文件
│   ├── multiLayerSolver.h          # 多层求解器头文件
│   ├── inverseAlgorithm.h          # 反演算法头文件
│   └── heatTransferSolverBase.h    # 基类头文件
└── build/                      # 构建输出目录 (自动生成)
    ├── bin/                    # 可执行程序
    ├── lib/                    # 库文件
    └── CMakeFiles/             # CMake临时文件
```

## 故障排除

### 问题1: 找不到CMake
**解决方案**: 在Ubuntu中安装CMake
```bash
sudo apt install cmake
```

### 问题2: 找不到Eigen3
**解决方案**: 安装Eigen3库
```bash
sudo apt install libeigen3-dev
```

### 问题3: 找不到Ceres Solver
**解决方案**: 安装Ceres库
```bash
sudo apt install libceres-dev
```

### 问题4: 编译错误 - 未定义的引用
**解决方案**: 确保所有依赖库都已正确安装，并且CMake能找到它们

### 问题5: 权限被拒绝 (build.sh)
**解决方案**: 添加执行权限
```bash
chmod +x build.sh
```

## 编译选项

### 优化编译
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -march=native" ..
```

### 调试编译
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

### 并行编译
```bash
cmake --build . --parallel $(nproc)
```

## 清理构建

```bash
# 删除构建目录
rm -rf build

# 或在build目录中
cd build
cmake --build . --target clean
```

## 支持的编译器

- **gcc** - 推荐
- **Clang** - 需要安装clang包
  ```bash
  sudo apt install clang
  ```

## 许可证


## 联系方式

