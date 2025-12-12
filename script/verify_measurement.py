#!/usr/bin/env python3
"""验证测量数据是否用真实值生成"""

import subprocess
import shutil
import numpy as np

print("="*70)
print("验证测量数据的生成")
print("="*70)

# 1. 备份当前measurement.txt
shutil.copy('measurement.txt', 'measurement_original.txt')

# 2. 使用真实值重新生成测量数据
print("\n步骤1: 使用真实热导率值生成新的测量数据...")

# 修改配置为mode=0（正演计算）
with open('config/single_layer_config.txt', 'r') as f:
    config = f.read()

config_forward = config.replace('SOLVE_MODE 1', 'SOLVE_MODE 0')

# 使用真实值
lines = config_forward.split('\n')
new_lines = []
in_material = False

for line in lines:
    if 'MATERIAL_DATA_START' in line:
        in_material = True
        new_lines.append(line)
    elif 'MATERIAL_DATA_END' in line:
        in_material = False
        new_lines.append(line)
    elif in_material and not line.startswith('#') and '\t' in line:
        parts = line.split('\t')
        if len(parts) >= 3:
            temp = parts[0]
            cp = parts[2]
            true_k = {'25': 50, '50': 52, '100': 54, '150': 56, '200': 58,
                     '250': 59, '350': 60, '400': 61, '450': 62, '500': 63, '550': 63}
            k = true_k.get(temp, 50)
            new_lines.append(f"{temp}\t{k}\t{cp}")
        else:
            new_lines.append(line)
    else:
        new_lines.append(line)

with open('config/single_layer_config_forward.txt', 'w') as f:
    f.write('\n'.join(new_lines))

# 运行正演计算
result = subprocess.run(['build/bin/singleLayerInverse', 'config/single_layer_config_forward.txt'],
                       capture_output=True, text=True, cwd='/home/dw/Opt')

# 重命名输出
shutil.move('temperature_output.txt', 'measurement_true.txt')

print("新测量数据已生成: measurement_true.txt")

# 3. 对比两个测量文件
print("\n步骤2: 对比原始测量数据和真实值生成的数据...")

data_orig = np.loadtxt('measurement_original.txt', skiprows=2)
data_true = np.loadtxt('measurement_true.txt', skiprows=2)

print(f"\n原始数据形状: {data_orig.shape}")
print(f"真实数据形状: {data_true.shape}")

# 对比T1列（索引2）
T1_orig = data_orig[:, 2]
T1_true = data_true[:, 2]

diff = np.abs(T1_orig - T1_true)
max_diff = np.max(diff)
mean_diff = np.mean(diff)
rel_diff = mean_diff / np.mean(T1_true) * 100

print(f"\nT1测点温度对比:")
print(f"  最大差异: {max_diff:.4f} °C")
print(f"  平均差异: {mean_diff:.4f} °C")
print(f"  相对差异: {rel_diff:.4f}%")

if max_diff < 0.1:
    print("\n✓ 测量数据确实是用真实值生成的")
else:
    print("\n✗ 测量数据NOT是用真实值生成的！")
    print("  这可能是误差的主要来源")
    
    # 显示几个时间点的对比
    print("\n时间点对比 (前10个):")
    print("Time\tT1_orig\tT1_true\tDiff")
    for i in range(min(10, len(T1_orig))):
        print(f"{data_orig[i,0]:.1f}\t{T1_orig[i]:.2f}\t{T1_true[i]:.2f}\t{diff[i]:.2f}")

# 恢复原文件
shutil.copy('measurement_original.txt', 'measurement.txt')
