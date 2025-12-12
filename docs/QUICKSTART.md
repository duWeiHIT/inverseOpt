# 快速开始指南 - WSL2 Ubuntu

## 5分钟快速构建

### 第1步: 打开WSL2 Ubuntu终端

### 第2步: 安装依赖 (首次运行)

```bash
sudo apt update
sudo apt install -y gcc g++ cmake make libeigen3-dev libceres-dev git
```

### 第3步: 进入项目目录

```bash
cd /path/to/Opt
```

### 第4步: 运行构建脚本

```bash
chmod +x build.sh
./build.sh
```

### 第5步: 运行程序

```bash
cd build/bin
./singleLayerInverse
```

## 常用命令

| 命令 | 说明 |
|------|------|
| `./build.sh` | 完整构建 |
| `cd build && cmake --build . --parallel $(nproc)` | 快速重新编译 |
| `rm -rf build` | 清理构建 |
| `./build/bin/singleLayerInverse` | 运行单层反演 |
| `./build/bin/multiLayerInverse` | 运行多层反演 |
| `./cmake_config.sh -d` | Debug构建 |
| `./check_dependencies.sh` | 检查依赖 |

## 环境变量设置 (可选)

如果CMake找不到库，可以手动指定:

```bash
# 如果需要指定Eigen3路径
export CMAKE_PREFIX_PATH="/usr/include/eigen3:$CMAKE_PREFIX_PATH"

# 或者设置PKG_CONFIG_PATH
export PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig:$PKG_CONFIG_PATH"
```

## 验证安装

```bash
# 检查CMake
cmake --version

# 检查编译器
gcc --version
g++ --version

# 检查Eigen3
pkg-config --modversion eigen3

# 检查Ceres
pkg-config --modversion ceres

# 一键检查所有依赖
./check_dependencies.sh
```

## WSL2特定提示

### 性能优化
```bash
# 检查CPU核数
nproc

# 使用所有核心编译
make -j$(nproc)

# 或限制并行数（如果内存不足）
make -j2
```

### 文件系统提示
```bash
# 如果项目在Windows分区，建议复制到WSL文件系统
cp -r /mnt/c/path/to/project ~/Opt
cd ~/Opt

# 修复权限问题
chmod +x *.sh
```

## 下一步

- 查看 `README.md` 了解详细信息
- 查看源代码了解程序功能
- 根据需要修改CMakeLists.txt
