#include "inverseAlgorithm.h"

int main() 
{
    try 
    { 
        std::cout<<"创建三层材料层"<<std::endl;
        std::vector<MaterialLayer> layers;
        for(unsigned int i = 0; i < 3; i++)
        {
            std::string material_name;
            double length;
            std::cout<<"输入样件"<<i<<"名字：";
            std::cin>>material_name;
            std::cout<<"输入样件: ["<<material_name<<"] 厚度(米)：";
            std::cin>>length;

            int N;
            std::cout<<"输入样件 ["<<material_name<<"] 空间离散点数目(外节点法): ";
            std::cin>>N;

            double rho;
            std::cout<<"输入样件 ["<<material_name<<"] 密度(kg/m^3): ";
            std::cin>>rho;
            
            layers.push_back(MaterialLayer(material_name, length, N, rho));

            std::cout<<"输入样件 ["<<material_name<<"] 热导率和比热(分段线性)"<<std::endl;
            
            int T_num;
            std::cout<<"输入温度点数目：";
            std::cin>>T_num;
            std::vector<double> Temperature(T_num);
            std::vector<double> kappa(T_num);
            std::vector<double> cp(T_num);

            std::cout<<"输入数据格式\n温度(摄氏度)\t热导率(W/m/K)\t比热(J/kg/K)"<<std::endl;
            for(int j = 0; j < T_num; j++)
            {
                std::cin>>Temperature[j]>>kappa[j]>>cp[j];
            }
            layers[i].thermal_conductivity.setData(Temperature, kappa);
            layers[i].specific_heat.setData(Temperature, cp);
        }
        
        BoundaryCondition left_bc;
        BoundaryCondition right_bc;
        
        std::cout<<"边界条件：\n"<<"0  定壁温\n"<<"1 时变定壁温\n"
        <<"2  定热流\n"<<"3  对流换热\n"<<"4  对流换热+辐射"<<std::endl;
        int left_type, right_type;

        while(1)
        {
            std::cout<<"输入左侧边界条件：";
            std::cin>>left_type;
            if(left_type == 0)
            {
                left_bc.type = BoundaryType::FIXED_TEMPERATURE;
                std::cout<<"输入左侧边界温度：";
                std::cin>>left_bc.value1;
                break;
            }
            else if(left_type == 1)
            {
                left_bc.type = BoundaryType::TIME_VARY_FIXED_TEMPERATURE;
                break;
            }
            else if(left_type == 2)
            {
                left_bc.type = BoundaryType::FIXED_HEAT_FLUX;
                std::cout<<"输入左侧边界热流：";
                std::cin>>left_bc.value1;
                break;
            }
            else if(left_type == 3)
            {
                left_bc.type = BoundaryType::CONVECTION;
                std::cout<<"输入左侧边界环境温度：";
                std::cin>>left_bc.value1;
                std::cout<<"输入左侧边界换热系数：";
                std::cin>>left_bc.value2;
                break;
            }
            else if(left_type == 4)
            {
                left_bc.type = BoundaryType::CONVECTION_RADIATION;
                std::cout<<"输入左侧边界环境温度：";
                std::cin>>left_bc.value1;
                std::cout<<"输入左侧边界换热系数：";
                std::cin>>left_bc.value2;
                std::cout<<"输入左侧边界发射率：";
                std::cin>>left_bc.value3;
                break;
            }
            else
            {
                std::cout<<"输入错误"<<std::endl;
                continue;
            }
        }

        while(1)
        {
            std::cout<<"输入右侧边界条件：";
            std::cin>>right_type;
            if(right_type == 0)
            {
                right_bc.type = BoundaryType::FIXED_TEMPERATURE;
                std::cout<<"输入右侧边界温度：";
                std::cin>>right_bc.value1;
                break;
            }
            else if(right_type == 1)
            {
                right_bc.type = BoundaryType::TIME_VARY_FIXED_TEMPERATURE;
                break;
            }
            else if(right_type == 2)
            {
                right_bc.type = BoundaryType::FIXED_HEAT_FLUX;
                std::cout<<"输入右侧边界热流：";
                std::cin>>right_bc.value1;
                break;
            }
            else if(right_type == 3)
            {
                right_bc.type = BoundaryType::CONVECTION;
                std::cout<<"输入右侧边界环境温度：";
                std::cin>>right_bc.value1;
                std::cout<<"输入右侧边界换热系数：";
                std::cin>>right_bc.value2;
                break;
            }
            else if(right_type == 4)
            {
                right_bc.type = BoundaryType::CONVECTION_RADIATION;
                std::cout<<"输入右侧边界环境温度：";
                std::cin>>right_bc.value1;
                std::cout<<"输入右侧边界换热系数：";
                std::cin>>right_bc.value2;
                std::cout<<"输入右侧边界发射率：";
                std::cin>>right_bc.value3;
                break;
            }
            else
            {
                std::cout<<"输入错误"<<std::endl;
                continue;
            }
        }

        double time_step;
        std::cout<<"输入时间步长, dt = ";
        std::cin>>time_step;
        
        std::cout<<"创建求解器"<<std::endl;
        multiLayerSolver1D solver(layers, time_step, left_bc, right_bc);
        
        // 打印层信息
        solver.printLayerInfo();
        
        // 设置初始温度
        double T_ini;
        std::cout<<"初始化....\n"<<"输入初始温度（摄氏度）：";
        std::cin>>T_ini;
        solver.initialise(T_ini);
        std::cout<<"初始化完成"<<std::endl;

        std::cout<<"导热计算"<<std::endl;
        std::cout<<"输入导热计算时间, t = ";
        double duration;
        std::cin>>duration;
        int time_step_num = std::round(duration/time_step);
        
        std::cout<<"输入监测点数目：";
        unsigned int measure_num;
        std::cin>>measure_num;
        std::vector<int> measure_index(measure_num);
        for(unsigned int i = 0; i < measure_num; i++)
        {
            std::cout<<"输入监测点"<<i<<"的位置, x = "<<std::endl;
            double measure_pos;
            std::cin>>measure_pos;
            const auto& node_positions = solver.getNodePositions();
            double min_distance = std::abs(node_positions[0] - measure_pos);
            for(int j = 0; j < solver.getTotalNodes(); j++)
            {
                double distance = std::abs(node_positions[j] - measure_pos);
                if (distance < min_distance) 
                {
                    min_distance = distance;
                    measure_index[i] = j;
                }
            }
        }
        auto& position = solver.getNodePositions();
        std::ofstream fout("result.txt");
        //写入文件头
        fout<<"time";
        for(unsigned int i = 0; i < measure_num; i++)
        {
            fout<<"\tT"<<i<<" pos ="<<position[measure_index[i]];
        }

        fout<<std::endl;
 
        for(int i = 0; i < time_step_num; i++)
        {
            solver.solve();
            auto& temperature = solver.getTemperature();
            fout<<i*time_step;
            for(unsigned int j = 0; j < measure_num; j++)
            {
                fout<<"\t"<<temperature[measure_index[j]];
            }

            fout<<std::endl;
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}