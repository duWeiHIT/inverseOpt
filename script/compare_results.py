#!/usr/bin/env python3
"""比较反演结果与真实值"""

import sys

def read_results(filename):
    """读取结果文件"""
    temps = []
    k_values = []
    
    with open(filename, 'r') as f:
        in_k_section = False
        for line in f:
            line = line.strip()
            if line.startswith('# 热导率反演结果'):
                in_k_section = True
                continue
            if line.startswith('# 比热'):
                break
            if in_k_section and line and not line.startswith('#'):
                parts = line.split()
                if len(parts) == 2:
                    temps.append(float(parts[0]))
                    k_values.append(float(parts[1]))
    
    return temps, k_values

def main():
    # 读取真实值和反演结果
    temps_true, k_true = read_results('result_true.txt')
    temps_inv, k_inv = read_results('inversion_result.txt')
    
    print("=" * 70)
    print("热导率反演结果对比")
    print("=" * 70)
    print(f"{'温度(°C)':<12} {'真实值':<15} {'反演值':<15} {'误差':<15} {'相对误差(%)':<15}")
    print("-" * 70)
    
    total_error = 0
    max_error = 0
    
    for i in range(len(temps_true)):
        error = abs(k_inv[i] - k_true[i])
        rel_error = error / k_true[i] * 100
        total_error += error
        max_error = max(max_error, error)
        
        print(f"{temps_true[i]:<12.0f} {k_true[i]:<15.2f} {k_inv[i]:<15.4f} {error:<15.4f} {rel_error:<15.2f}")
    
    print("-" * 70)
    print(f"平均绝对误差: {total_error/len(temps_true):.4f} W/m/K")
    print(f"最大绝对误差: {max_error:.4f} W/m/K")
    print(f"平均相对误差: {total_error/len(temps_true)/sum(k_true)*len(temps_true)*100:.2f}%")
    print("=" * 70)

if __name__ == '__main__':
    main()
