#!/bin/bash
# 测试不同的正则化参数

echo "测试不同正则化参数的效果"
echo "======================================"

lambdas=(0 0.0005 0.001 0.003 0.005 0.008 0.01 0.03 0.1)

for lambda in "${lambdas[@]}"; do
    echo ""
    echo "测试 λ = $lambda"
    echo "--------------------------------------"
    
    # 修改配置文件
    sed -i "s/^REGULARIZATION_LAMBDA.*/REGULARIZATION_LAMBDA $lambda/" config/single_layer_config.txt
    
    # 运行反演
    build/bin/singleLayerInverse > /dev/null 2>&1
    
    # 显示结果
    echo "反演结果:"
    python3 compare_results.py | tail -5
    
    # 保存结果
    cp inversion_result.txt "inversion_result_lambda_${lambda}.txt"
done

echo ""
echo "======================================"
echo "所有测试完成，结果已保存"
