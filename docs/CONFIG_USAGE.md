# 配置文件使用说明

## 概述
现在程序支持从配置文件读取参数，避免了繁琐的交互式输入。

## 文件结构
```
config/
├── single_layer_config.txt     # 单层材料配置
├── multi_layer_config.txt      # 多层材料配置
└── measurement.txt             # 示例测量数据
```

## 使用方法

### 单层材料反演
```bash
# 使用默认配置文件
./build/bin/singleLayerInverse

# 使用自定义配置文件
./build/bin/singleLayerInverse path/to/your/config.txt
```

### 多层材料反演
```bash
# 使用默认配置文件
./build/bin/multiLayerInverse

# 使用自定义配置文件
./build/bin/multiLayerInverse path/to/your/config.txt
```

## 配置文件格式

### 单层材料配置 (single_layer_config.txt)
```
# 几何参数
LENGTH 0.01                    # 样件长度(m)
NODES 50                       # 离散单元数目

# 边界条件类型: 0=定壁温 1=时变定壁温 2=定热流 3=对流换热 4=对流换热+辐射
LEFT_BC_TYPE 0
LEFT_BC_VALUE1 100.0           # 参数1(温度/热流/环境温度)
LEFT_BC_VALUE2 0.0             # 参数2(换热系数)
LEFT_BC_VALUE3 0.0             # 参数3(发射率)

RIGHT_BC_TYPE 3
RIGHT_BC_VALUE1 25.0
RIGHT_BC_VALUE2 10.0
RIGHT_BC_VALUE3 0.0

# 材料属性
DENSITY 7800.0                 # 密度(kg/m^3)

# 材料物性数据
MATERIAL_DATA_START
25.0  50.0  500.0             # 温度(°C) 热导率(W/m/K) 比热(J/kg/K)
100.0 45.0  520.0
200.0 40.0  550.0
MATERIAL_DATA_END

# 时间参数
TIME_STEP 0.1                  # 时间步长(s)
INITIAL_TEMP 25.0              # 初始温度(°C)

# 求解参数
SOLVE_MODE 1                   # 0=导热计算 1=物性反演
DURATION 10.0                  # 导热计算时间(s)

# 反演参数
MEASUREMENT_FILE "measurement.txt"
MEASURE_POSITION 0.005         # 测点位置(m)
ALGORITHM_INDEX 0              # 算法索引
USE_DEFAULT_CONFIG true        # 使用默认算法配置
```

### 多层材料配置 (multi_layer_config.txt)
```
# 层数
LAYER_COUNT 3

# 各层参数
LAYER1_NAME "Material1"
LAYER1_THICKNESS 0.003
LAYER1_NODES 20
LAYER1_DENSITY 7800.0
LAYER1_TEMP_POINTS 3
LAYER1_DATA_START
25.0  50.0  500.0
100.0 45.0  520.0
200.0 40.0  550.0
LAYER1_DATA_END

# ... 其他层类似定义

# 反演参数
SAMPLE_INDEX 0                 # 待反演样件编号
```

## 边界条件类型说明
- 0: 定壁温 - 需要VALUE1(温度)
- 1: 时变定壁温 - 程序内部处理
- 2: 定热流 - 需要VALUE1(热流)
- 3: 对流换热 - 需要VALUE1(环境温度), VALUE2(换热系数)
- 4: 对流换热+辐射 - 需要VALUE1(环境温度), VALUE2(换热系数), VALUE3(发射率)

## 注意事项
1. 配置文件中的注释以#开头
2. 参数名称区分大小写
3. 字符串值需要用双引号包围
4. 布尔值使用true/false
5. 测量数据文件路径相对于程序运行目录

## 新增功能

### 温度输出功能 (mode=0)
当求解模式设为0（导热计算）时，程序会自动输出监测点的温度数据到文件。

#### 配置参数
```
OUTPUT_FILE "temperature_output.txt"  # 输出文件名
MONITOR_POSITION 0.005               # 监测点位置(m)
OUTPUT_INTERVAL 1                    # 输出间隔
```

#### 输出文件格式
```
# 温度测量数据文件
# 监测点位置: 0.005 m
# 格式: 时间(s) 温度(°C)
0.0  25.0
0.1  28.5
0.2  32.1
...
```

#### 使用场景
- **程序验证**: 生成模拟数据用于反演程序验证
- **参数研究**: 分析不同参数对温度响应的影响
- **结果对比**: 与实验数据或其他模型结果对比

## 优势
- 避免重复输入参数
- 便于批量计算和参数研究
- 配置可重用和版本控制
- 减少输入错误
- 自动温度数据输出