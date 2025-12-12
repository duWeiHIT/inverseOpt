#!/usr/bin/env python3
"""测试使用真实值作为初始猜测"""

import subprocess
import shutil

# 备份原配置
shutil.copy('config/single_layer_config.txt', 'config/single_layer_config.txt.bak')

# 修改配置文件，使用真实值作为初始猜测
with open('config/single_layer_config.txt', 'r') as f:
    lines = f.readlines()

with open('config/single_layer_config.txt', 'w') as f:
    in_material_section = False
    for line in lines:
        if 'MATERIAL_DATA_START' in line:
            in_material_section = True
            f.write(line)
        elif 'MATERIAL_DATA_END' in line:
            in_material_section = False
            f.write(line)
        elif in_material_section and not line.startswith('#'):
            # 使用真实值
            parts = line.split()
            if len(parts) == 3:
                temp = float(parts[0])
                cp = parts[2]
                # 真实热导率值
                true_k = {25: 50, 50: 52, 100: 54, 150: 56, 200: 58, 
                         250: 59, 350: 60, 400: 61, 450: 62, 500: 63, 550: 63}
                k = true_k.get(temp, 50)
                f.write(f"{temp}\t{k}\t{cp}\n")
            else:
                f.write(line)
        else:
            f.write(line)

print("已修改配置文件，使用真实值作为初始猜测")
print("运行反演...")

# 运行反演
result = subprocess.run(['build/bin/singleLayerInverse'], 
                       capture_output=True, text=True)

print("\n" + "="*70)
print("使用真实值作为初始猜测的反演结果:")
print("="*70)

# 显示结果
with open('inversion_result.txt', 'r') as f:
    print(f.read())

# 恢复原配置
shutil.copy('config/single_layer_config.txt.bak', 'config/single_layer_config.txt')
print("\n已恢复原配置文件")
