

#include "inverseAlgorithm.h"
#include "configReader.h"

// 从配置文件读取材料数据
std::vector<std::vector<double>> readMaterialDataFromConfig(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;
    std::string line;
    bool reading_data = false;
    
    while (std::getline(file, line)) {
        if (line == "MATERIAL_DATA_START") {
            reading_data = true;
            continue;
        }
        if (line == "MATERIAL_DATA_END") {
            break;
        }
        if (reading_data && !line.empty() && line[0] != '#') {
            std::istringstream iss(line);
            std::vector<double> row;
            double value;
            while (iss >> value) {
                row.push_back(value);
            }
            if (row.size() == 3) {
                data.push_back(row);
            }
        }
    }
    return data;
}

// 使用示例
int main(int argc, char* argv[]) 
{
    try 
    {
        std::string config_file = "config/single_layer_config.txt";
        if (argc > 1) {
            config_file = argv[1];
        }
        
        ConfigReader config;
        if (!config.loadConfig(config_file)) {
            std::cerr << "配置文件加载失败，使用交互式输入" << std::endl;
            // 保留原有交互式输入代码作为备选
        }
        
        double L = config.getDouble("LENGTH");
        int nodes = config.getInt("NODES");
        
        std::cout << "从配置文件读取: L = " << L << ", nodes = " << nodes << std::endl;

        BoundaryCondition left_bc;
        BoundaryCondition right_bc;
        
        // 从配置文件读取边界条件
        int left_type = config.getInt("LEFT_BC_TYPE");
        left_bc.value1 = config.getDouble("LEFT_BC_VALUE1");
        left_bc.value2 = config.getDouble("LEFT_BC_VALUE2");
        left_bc.value3 = config.getDouble("LEFT_BC_VALUE3");
        
        switch(left_type) {
            case 0: left_bc.type = BoundaryType::FIXED_TEMPERATURE; break;
            case 1: left_bc.type = BoundaryType::TIME_VARY_FIXED_TEMPERATURE; break;
            case 2: left_bc.type = BoundaryType::FIXED_HEAT_FLUX; break;
            case 3: left_bc.type = BoundaryType::CONVECTION; break;
            case 4: left_bc.type = BoundaryType::CONVECTION_RADIATION; break;
        }
        
        int right_type = config.getInt("RIGHT_BC_TYPE");
        right_bc.value1 = config.getDouble("RIGHT_BC_VALUE1");
        right_bc.value2 = config.getDouble("RIGHT_BC_VALUE2");
        right_bc.value3 = config.getDouble("RIGHT_BC_VALUE3");
        
        switch(right_type) {
            case 0: right_bc.type = BoundaryType::FIXED_TEMPERATURE; break;
            case 1: right_bc.type = BoundaryType::TIME_VARY_FIXED_TEMPERATURE; break;
            case 2: right_bc.type = BoundaryType::FIXED_HEAT_FLUX; break;
            case 3: right_bc.type = BoundaryType::CONVECTION; break;
            case 4: right_bc.type = BoundaryType::CONVECTION_RADIATION; break;
        }
        
        std::cout << "边界条件已从配置文件读取" << std::endl;

        double rho = config.getDouble("DENSITY");
        
        std::cout << "从配置文件读取材料属性: 密度 = " << rho << " kg/m^3" << std::endl;
        
        // 读取材料数据
        auto material_data = readMaterialDataFromConfig(config_file);
        int T_num = material_data.size();
        
        std::vector<double> Temperature(T_num);
        std::vector<double> kappa(T_num);
        std::vector<double> cp(T_num);
        
        for(int j = 0; j < T_num; j++) 
        {
            Temperature[j] = material_data[j][0];
            kappa[j] = material_data[j][1];
            cp[j] = material_data[j][2];
        }
        
        std::cout << "材料数据点数: " << T_num << std::endl;

        MaterialLayer layer("single_layer", L, nodes, rho);
        layer.thermal_conductivity.setData(Temperature, kappa);
        layer.specific_heat.setData(Temperature, cp);

        MaterialProperty& kappa_prop = layer.thermal_conductivity;
        MaterialProperty& cp_prop = layer.specific_heat;

        double time_step = config.getDouble("TIME_STEP");
        double T_ini = config.getDouble("INITIAL_TEMP");
        int mode = config.getInt("SOLVE_MODE");
        
        std::cout << "时间步长: " << time_step << " s" << std::endl;
        std::cout << "初始温度: " << T_ini << " °C" << std::endl;
        std::cout << "求解模式: " << (mode == 0 ? "导热计算" : "物性反演") << std::endl;

        std::cout<<"初始化求解器"<<std::endl;
        singleLayerSolver1D solver(layer, time_step, left_bc, right_bc);
        solver.initialise(T_ini);

        std::cout<<"初始化完成"<<std::endl;
        if(mode == 0)
        {
            std::cout<<"导热计算"<<std::endl;
            double duration = config.getDouble("DURATION");
            std::cout << "计算时间: " << duration << " s" << std::endl;
            
            // 读取输出参数
            std::string output_file = config.getString("OUTPUT_FILE");
            int monitor_count = config.getInt("MONITOR_COUNT");
            int output_interval = config.getInt("OUTPUT_INTERVAL");
            
            // 读取多个监测点位置
            std::vector<double> monitor_positions(monitor_count);
            std::vector<int> monitor_indices(monitor_count);
            const double dx = L/(nodes-1);
            
            for(int i = 0; i < monitor_count; i++)
            {
                std::string pos_key = "MONITOR_POSITION" + std::to_string(i + 1);
                monitor_positions[i] = config.getDouble(pos_key);
                
                // 找到最近的节点
                double min_distance = std::abs(0 - monitor_positions[i]);
                monitor_indices[i] = 0;
                for(int j = 0; j < nodes; j++)
                {
                    double distance = std::abs(j*dx - monitor_positions[i]);
                    if(distance < min_distance)
                    {
                        min_distance = distance;
                        monitor_indices[i] = j;
                    }
                }
                
                std::cout << "监测点" << i+1 << ": " << monitor_positions[i] << " m (节点" << monitor_indices[i] << ", 实际位置: " << monitor_indices[i]*dx << " m)" << std::endl;
            }
            std::cout << "输出文件: " << output_file << std::endl;
            
            // 创建输出文件
            std::ofstream fout(output_file);
            if (!fout.is_open()) {
                std::cerr << "无法创建输出文件: " << output_file << std::endl;
                return 1;
            }
            
            // 写入文件头
            fout << "# 温度测量数据文件" << std::endl;
            fout << "time";
            for(int i = 0; i < monitor_count; i++)
            {
                fout << "\tT" << i;
            }
            fout << std::endl;
            
            int time_step_num = std::round(duration/time_step);
            
            // 输出初始温度
            const auto& temperature = solver.getTemperature();
            fout << "0.0";
            for(int i = 0; i < monitor_count; i++)
            {
                fout << "\t" << temperature[monitor_indices[i]];
            }
            fout << std::endl;
            
            for(int i = 0; i < time_step_num; i++)
            {
                solver.solve();
                
                // 按间隔输出
                if((i + 1) % output_interval == 0)
                {
                    const auto& temp = solver.getTemperature();
                    double current_time = (i + 1) * time_step;
                    fout << current_time;
                    for(int j = 0; j < monitor_count; j++)
                    {
                        fout << "\t" << temp[monitor_indices[j]];
                    }
                    fout << std::endl;
                }
            }
            
            fout.close();
            std::cout << "温度数据已输出到: " << output_file << std::endl;
        }
        else if(mode == 1)
        {
            std::cout<<"物性反演"<<std::endl;
            std::string filename = config.getString("MEASUREMENT_FILE");
            int inversion_mode = config.getInt("INVERSION_MODE");
            if(inversion_mode == 0)
            {
                std::cout << "反演模式: 同时反演热导率和比热" << std::endl;
            }
            else if(inversion_mode == 1)
            {
                std::cout << "反演模式: 仅反演热导率" << std::endl;
            }
            else
            {
                throw std::runtime_error("无效的反演模式");
            }

            int measure_point_count = config.getInt("MEASURE_POINT_COUNT");
            
            std::cout << "温度测量文件: " << filename << std::endl;
            std::cout << "测点数量: " << measure_point_count << std::endl;
            
            std::vector<int> measure_positions;
            std::vector<int> measure_indices;
            const double dx = L/(nodes-1);
            
            for(int p = 1; p <= measure_point_count; p++) 
            {
                std::string pos_key = "MEASURE_POSITION" + std::to_string(p);
                std::string idx_key = "MEASURE_POINT_INDEX" + std::to_string(p);
                
                double measure_pos = config.getDouble(pos_key);
                int measure_point_index = config.getInt(idx_key);
                
                int measure_node = 0;
                double min_distance = std::abs(0 - measure_pos);
                for(int i = 0; i < nodes; i++) 
                {
                    double distance = std::abs(i*dx - measure_pos);
                    if(distance < min_distance) 
                    {
                        min_distance = distance;
                        measure_node = i;
                    }
                }
                
                measure_positions.push_back(measure_node);
                measure_indices.push_back(measure_point_index);

                std::cout << "测点" << p << ": 位置=" << measure_pos << "m, 编号=" << measure_point_index << std::endl;
            }

            std::cout<<"读取温度测量文件: "<<filename<<std::endl;
            solver.loadMeasurementData(filename, measure_positions, measure_indices);
            std::cout<<"total step = "<<solver.getTotalSteps();
            std::cout<<"测量文件读取完成"<<std::endl;
            std::cout<<"创建反演目标函数"<<std::endl;
            ThermalPropertyObjective objective(solver, kappa_prop, cp_prop, inversion_mode);
            
            double k_lower_val = config.getDouble("K_LOWER_BOUND");
            double k_upper_val = config.getDouble("K_UPPER_BOUND");
            double cp_lower_val = config.getDouble("CP_LOWER_BOUND");
            double cp_upper_val = config.getDouble("CP_UPPER_BOUND");
            
            objective.setThermalConductivityBounds(k_lower_val, k_upper_val);
            objective.setSpecificHeatBounds(cp_lower_val, cp_upper_val);
            
            double reg_lambda = config.getDouble("REGULARIZATION_LAMBDA");
            objective.setRegularizationParameter(reg_lambda);
            
            std::cout << "参数边界设置: 热导率[" << k_lower_val << "-" << k_upper_val 
                      << "] W/m/K, 比热[" << cp_lower_val << "-" << cp_upper_val << "] J/kg/K" << std::endl;
            std::cout << "正则化参数: " << reg_lambda << std::endl;
            
            std::vector<double> initial_guess(objective.getParameterCount());
            for(size_t i = 0; i < kappa_prop.getValues().size(); i++)
            {
                initial_guess[i] = kappa_prop.getValues()[i];
            }
            
            if(inversion_mode == 0) 
            {
                size_t k_count = kappa_prop.getValues().size();
                for(size_t i = 0; i < cp_prop.getValues().size(); i++)
                {
                    initial_guess[i + k_count] = cp_prop.getValues()[i];
                }
            }

            std::cout << "\n初始猜测值:" << std::endl;
            for(size_t i = 0; i < initial_guess.size(); i++) 
            {
                std::cout << "  param[" << i << "] = " << initial_guess[i] << std::endl;
            }

            //归一化
            std::cout << "\n归一化猜测值:" << std::endl;
            objective.normalizeParameters(initial_guess);

            // 创建算法工厂
            InversionAlgorithmFactory factory;
                    
            // 获取可用算法
            auto algorithms = factory.getAvailableAlgorithms();
            int index = config.getInt("ALGORITHM_INDEX");
            
            std::cout << "可用反演算法: ";
            for (size_t i = 0; i < algorithms.size(); i++) 
            {
                std::cout << i << " -- " << algorithms[i] << "\n";
            }
            std::cout << "选择的算法索引: " << index << std::endl;
            
            if(index < 0 || static_cast<size_t>(index) >= algorithms.size()) {
                std::cerr << "算法索引超出范围，使用默认算法 1" << std::endl;
                index = 1;
            }
            
            std::string algorithm_name = algorithms[index];
            
            std::cout << "\n=== 使用算法: " << algorithm_name << " ===" << std::endl;
            bool default_config = config.getBool("USE_DEFAULT_CONFIG");
            auto algorithm = factory.createAlgorithm(algorithm_name);
            
            std::cout<<"默认配置: "<<std::endl;
            algorithm->outputConfig();
            
            if(default_config)
            {
                std::cout<<"使用默认配置"<<std::endl;
            }
            else
            {
                std::cout<<"使用自定义配置"<<std::endl;
                algorithm->setConfig();
            }
                
            auto result = algorithm->optimize(&objective, initial_guess);
            
            std::cout << "算法: " << algorithm->getName() << std::endl;
            std::cout << "最优值: " << result.optimal_value << std::endl;
            std::cout << "迭代次数: " << result.iterations << std::endl;
            std::cout << "是否收敛: " << (result.converged ? "是" : "否") << std::endl;
            std::cout << "消息: " << result.message << std::endl;
            objective.updateMaterialProperties(result.optimal_parameters);
            std::cout<< "优化参数：\n";
            for(size_t i = 0; i < result.optimal_parameters.size(); i++)
            {
                std::cout<<result.optimal_parameters[i]<<"\n";
            }
            std::cout<<std::endl;
            std::ofstream result_file("inversion_result.txt");
            
            if (result_file.is_open()) 
            {
                result_file << "# 热物性反演结果" << std::endl;
                result_file << "# 算法: " << algorithm->getName() << std::endl;
                result_file << "# 最优值: " << result.optimal_value << std::endl;
                result_file << "# 迭代次数: " << result.iterations << std::endl;
                result_file << "# 收敛状态: " << (result.converged ? "收敛" : "未收敛") << std::endl;
                result_file << "# 参数化模式: " << "分段线性" << std::endl;
                result_file << std::endl;
                

                result_file << "# 热导率反演结果" << std::endl;
                result_file << "# 温度(°C)\t热导率(W/m/K)" << std::endl;
                for (size_t i = 0; i < kappa_prop.getTemperatures().size(); i++)
                {
                    result_file << kappa_prop.getTemperatures()[i] << "\t" 
                                << kappa_prop.getValues()[i] << std::endl;
                }

                
                if(inversion_mode == 0) 
                {
                    result_file << std::endl << "# 比热反演结果" << std::endl;
                } 
                else 
                {
                    result_file << std::endl << "# 比热保持原值（未反演）" << std::endl;
                }

                result_file << "# 温度(°C)\t比热(J/kg/K)" << std::endl;
                for (size_t i = 0; i < cp_prop.getTemperatures().size(); i++) 
                {
                    result_file << cp_prop.getTemperatures()[i] << "\t" 
                                << cp_prop.getValues()[i] << std::endl;
                }
                
                result_file.close();
                std::cout << "反演结果已保存到: inversion_result.txt" << std::endl;
            } 
            else 
            {
                std::cerr << "无法创建结果文件" << std::endl;
            }
        }
        else
        {
            throw std::runtime_error("无效的求解模式");
        }
    } 
    catch (const std::exception& e) 
    {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}