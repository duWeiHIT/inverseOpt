#ifndef HEAT_TRANSFER_SOLVER_BASE_H
#define HEAT_TRANSFER_SOLVER_BASE_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <ceres/ceres.h>

//Stefan-Boltamann常数
const double stefan_boltzmann = 5.670367e-8;

// 边界条件类型枚举
enum class BoundaryType 
{
    FIXED_TEMPERATURE,           // 第一类：定壁温
    TIME_VARY_FIXED_TEMPERATURE, // 第二类：时变定壁温
    FIXED_HEAT_FLUX,             // 第三类：定热流
    CONVECTION,                  // 第四类：对流换热
    CONVECTION_RADIATION         // 第五类：对流换热+辐射
};

// 边界条件结构体
struct BoundaryCondition 
{
    BoundaryType type;
    double value1;  // 温度或热流密度或流体温度
    double value2;  // 换热系数或辐射参数
    double value3;  // 辐射参数（发射率等)
};

// 材料属性类（支持分段线性和多项式）
class MaterialProperty 
{
private:
    std::vector<double> temperatures;
    std::vector<double> values; // 分段线性模式:热导率值; 多项式模式:系数
    std::string name;

public:
    MaterialProperty(const std::string& prop_name) 
    : name(prop_name) 
    {}
    
    void setData(const std::vector<double>& temps, const std::vector<double>& vals) 
    {
        if(temps.size() != vals.size())
        {
            throw std::invalid_argument("温度点和数值点数量不匹配");
        }
        temperatures = temps;
        values = vals;
    }

    void addDataPoint(double temp, double val) 
    {
        // 找到插入位置以保持温度有序
        auto it = std::lower_bound(temperatures.begin(), temperatures.end(), temp);
        size_t pos = it - temperatures.begin();
        
        temperatures.insert(it, temp);
        values.insert(values.begin() + pos, val);
    }

    // void autoSegment
    // (
    //     const std::vector<double>& temp_data, 
    //     const std::vector<double>& val_data, 
    //     int num_segments
    // ) 
    // {
    //     if (temp_data.size() != val_data.size() || temp_data.empty()) 
    //     {
    //         throw std::invalid_argument("输入数据无效");
    //     }
        
    //     // 找到温度范围
    //     double min_temp = *std::min_element(temp_data.begin(), temp_data.end());
    //     double max_temp = *std::max_element(temp_data.begin(), temp_data.end());
        
    //     // 清空现有数据
    //     temperatures.clear();
    //     values.clear();
        
    //     // 如果分段数大于数据点数，使用所有数据点
    //     if (static_cast<size_t>(num_segments) >= temp_data.size()) 
    //     {
    //         temperatures = temp_data;
    //         values = val_data;
    //         return;
    //     }
        
    //     // 均匀分段
    //     double segment_size = (max_temp - min_temp) / num_segments;
    //     for (int i = 0; i <= num_segments; ++i) 
    //     {
    //         double temp = min_temp + i * segment_size;
            
    //         // 找到该温度附近的点进行线性回归
    //         std::vector<std::pair<double, double>> nearby_points;
    //         for (size_t j = 0; j < temp_data.size(); ++j) {
    //             if (std::abs(temp_data[j] - temp) <= segment_size) {
    //                 nearby_points.emplace_back(temp_data[j], val_data[j]);
    //             }
    //         }
            
    //         if (!nearby_points.empty()) 
    //         {
    //             // 计算平均值作为该分段点的值
    //             double sum_val = 0.0;
    //             for (const auto& point : nearby_points) {
    //                 sum_val += point.second;
    //             }
    //             double avg_val = sum_val / nearby_points.size();
                
    //             temperatures.push_back(temp);
    //             values.push_back(avg_val);
    //         }
    //     }
    // }
    
    double evaluate(double temp) const 
    {
        if (temperatures.empty()) return 0.0;
        if (temp <= temperatures.front()) return values.front();
        if (temp >= temperatures.back()) return values.back();
        
        for (size_t i = 0; i < temperatures.size() - 1; i++) 
        {
            if (temp >= temperatures[i] && temp <= temperatures[i + 1]) 
            {
                double t1 = temperatures[i];
                double t2 = temperatures[i + 1];
                double v1 = values[i];
                double v2 = values[i + 1];

                if (std::abs(t2 - t1) < 1e-10) 
                {
                    return (v1 + v2) / 2.0;
                }

                return v1 + (v2 - v1) * (temp - t1) / (t2 - t1);
            }
        }
        return values.back();
    }
    
    void print() const 
    {

        std::cout << name << " 分段线性函数:" << std::endl;
        for (size_t i = 0; i < temperatures.size(); ++i) 
        {
            std::cout << "  T = " << temperatures[i] << "K -> " << values[i];
            if (i < temperatures.size() - 1) 
            {
                std::cout << " (斜率: " << (values[i+1] - values[i]) / (temperatures[i+1] - temperatures[i]) << ")";
            }
            std::cout << std::endl;
        }
    }

    size_t getNumSegments() const 
    {
        return temperatures.size();
    }

    const std::vector<double>& getTemperatures() const 
    {
        return temperatures;
    }

    const std::vector<double>& getValues() const 
    {
        return values;
    }

    std::vector<double>& getValues()
    {
        return values;
    }

    std::string getName() const { return name; }
};

// 自动微分兼容的材料属性类
template<typename T>
class MaterialPropertyAD 
{
private:
    std::vector<double> temperatures;
    std::vector<T> values;
    std::string name;

public:
    MaterialPropertyAD(const std::string& prop_name) 
    : name(prop_name) 
    {}
    
    void setData(const std::vector<double>& temps) 
    {
        temperatures = temps;
        values.resize(temps.size());
    }

    void setValue(const std::vector<T>& vals) 
    {
        if(temperatures.size() != vals.size())
        {
            throw std::invalid_argument("温度点和数值点数量不匹配");
        }

        values = vals;
    }

    void addDataPoint(double temp, double val) 
    {
        const T temp_vals = T(temp);
        auto it = std::lower_bound(temperatures.begin(), temperatures.end(), temp_vals);
        size_t pos = it - temperatures.begin();
        
        temperatures.insert(it, temp_vals);
        values.insert(values.begin() + pos, T(val));
    }

    T evaluate(T temp) const 
    {
        if (temperatures.empty()) return T(0.0);
        if (temp <= temperatures.front()) return values.front();
        if (temp >= temperatures.back()) return values.back();
        
        for (size_t i = 0; i < temperatures.size() - 1; i++) 
        {
            if (temp >= temperatures[i] && temp <= temperatures[i + 1]) 
            {
                const double t1 = temperatures[i];
                const double t2 = temperatures[i + 1];
                const T& v1 = values[i];
                const T& v2 = values[i + 1];

                if (std::abs(t2 - t1) < 1e-10) 
                {
                    return 0.5*(v1 + v2);
                }

                return v1 + (v2 - v1)*(temp - t1)/(t2 - t1);
            }
        }

        return values.back();
    }

    const std::vector<T>& getValues() const 
    {
        return values;
    }

    std::vector<T>& getValues() 
    {
        return values;
    }
};

// 材料层结构
struct MaterialLayer 
{
    double thickness;                       // 层厚度
    double density;                         // 密度
    int    nodes;                           // 离散点数目
    MaterialProperty thermal_conductivity;  // 热导率
    MaterialProperty specific_heat;         // 比热
    std::string name;                       // 材料名称

    MaterialLayer() 
    :
    thickness(0.01),
    density(1000),
    nodes(10),
    thermal_conductivity("k_default"),
    specific_heat("cp_default"),
    name("default")
    {}
    
    MaterialLayer(const std::string& layer_name, double thick, int n, double rho) 
    :
    thickness(thick),
    density(rho),
    nodes(n),
    thermal_conductivity("k_" + layer_name),
    specific_heat("cp_" + layer_name),
    name(layer_name)
    {}
    
    inline double getK(double temp) const 
    {
        return thermal_conductivity.evaluate(temp);
    }
    
    inline double getCp(double temp) const 
    {
        return specific_heat.evaluate(temp);
    }
};

// 自动微分兼容的材料层
template<typename T>
struct MaterialLayerAD 
{
    double thickness;
    double density;
    int nodes;
    MaterialPropertyAD<T> thermal_conductivity;
    MaterialPropertyAD<T> specific_heat;
    std::string name;

    MaterialLayerAD(const MaterialLayer& original_layer)
    :
    thickness(original_layer.thickness), 
    density(original_layer.density), 
    nodes(original_layer.nodes),
    thermal_conductivity("k_" + original_layer.name),
    specific_heat("cp_" + original_layer.name),
    name(original_layer.name)
    {
        // 复制热导率温度点
        thermal_conductivity.setData
        (
            original_layer.thermal_conductivity.getTemperatures()
        );
        
        // 复制比热温度点
        specific_heat.setData
        (
            original_layer.specific_heat.getTemperatures()
        );
    }

    inline T getK(T temp) const 
    {
        return thermal_conductivity.evaluate(temp);
    }
    
    inline T getCp(T temp) const 
    {
        return specific_heat.evaluate(temp);
    }
};

// 温度测量数据
class TemperatureMeasurement 
{
private:
    std::vector<double> times;
    std::vector<std::string> position_names;
    std::vector<std::vector<double>> temperatures;
    std::vector<int> measure_positions_; // 测点位置
    std::vector<int> measure_indices_;   // 测点索引

public:

bool loadFromFile(const std::string& filename) 
{
    std::ifstream file;

    try
    {
        file.open(filename);
        if (!file.is_open()) {
            throw std::runtime_error("无法打开文件: " + filename);
        }

        times.clear();
        position_names.clear();
        temperatures.clear();

        std::string line;
        size_t line_number = 0;
        size_t num_of_positions = 0;
        while (std::getline(file, line)) 
        {
            if (line.empty() || line[0] == '#') continue;
            if (line_number == 0) 
            {
                std::istringstream title_stream(line);
                std::string header;
                
                // 跳过第一个标题 "Times"
                if (!(title_stream >> header)) 
                {
                    throw std::runtime_error("文件格式错误：无法读取标题行");
                }
                
                // 读取后续的温度位置标题，如 T0, T1, T2...
                std::string pos_name;
                while (title_stream >> pos_name) 
                {
                    position_names.push_back(pos_name);
                }
                num_of_positions = position_names.size();
                
                if (num_of_positions == 0) 
                {
                    throw std::runtime_error("文件格式错误：未找到温度位置标题");
                }
                
                std::cout << "检测到 " << num_of_positions << " 个温度位置: ";
                for (const auto& name : position_names) 
                {
                    std::cout << name << " ";
                }
                std::cout << std::endl;
                
                line_number++;
                continue;
            }

            // 处理数据行
            std::istringstream data_stream(line);
            double current_time;
            // 读取时间值
            if (!(data_stream >> current_time)) 
            {
                throw std::runtime_error("数据格式错误在第" + std::to_string(line_number+1) + "行: " + line);
            }
            times.push_back(current_time);

            // 读取该时间点下所有位置的温度值
            std::vector<double> current_temps;
            double temp_value;
            for (size_t i = 0; i < num_of_positions; ++i) 
            {
                if(!(data_stream >> temp_value)) 
                {
                    throw std::runtime_error(
                                                "第" + std::to_string(line_number+1) 
                                                + "行数据不完整，期望" + 
                                               std::to_string(num_of_positions) 
                                               + "个温度值"
                                            );
                }
                current_temps.push_back(temp_value);
            }
            
            temperatures.push_back(current_temps);
            line_number++;
        }

        if (times.empty()) 
        {
            throw std::runtime_error("文件中没有有效数据");
        }

        std::cout << "成功加载 " << times.size() << " 个时间点的数据，每个时间点包含 " 
                << num_of_positions << " 个位置的温度值" << std::endl;
        return true;
    }
    catch(const std::ios_base::failure& e)
    {
        std::cerr << "【致命错误】文件I/O操作失败: " << e.what() << std::endl;
        return false;
    }
    catch (const std::exception& e) 
    {
        std::cerr << "【致命错误】数据处理失败: " << e.what() << std::endl;
        return false;
    }
}

    // 获取指定位置和时间的温度
    inline double getTemperature(int timeIndex, int position_index) const 
    {
        if(position_index >= 0 && position_index < static_cast<int>(position_names.size()) &&
           timeIndex >= 0 && timeIndex < static_cast<int>(temperatures.size()))
        {
            return temperatures[timeIndex][position_index];
        }
        return 0.0;
    }
    
    // 根据时间进行分段线性插值获取温度
    inline double getTemperatureByTime(double time, int position_index) const 
    {
        if (times.empty() || position_index < 0 || position_index >= static_cast<int>(position_names.size())) {
            return 0.0;
        }
        
        // 边界情况处理
        if (time <= times.front()) {
            return temperatures.front()[position_index];
        }
        if (time >= times.back()) {
            return temperatures.back()[position_index];
        }
        
        // 找到时间区间进行线性插值
        for (size_t i = 0; i < times.size() - 1; ++i) {
            if (time >= times[i] && time <= times[i + 1]) {
                double t1 = times[i];
                double t2 = times[i + 1];
                double T1 = temperatures[i][position_index];
                double T2 = temperatures[i + 1][position_index];
                
                // 线性插值
                if (std::abs(t2 - t1) < 1e-10) {
                    return (T1 + T2) / 2.0;
                }
                return T1 + (T2 - T1) * (time - t1) / (t2 - t1);
            }
        }
        
        return 0.0;
    }

    // 获取特定位置在所有时间点的温度
    std::vector<double> getTemperaturesAtPosition(int position_index) const 
    {
        std::vector<double> result;
        if(position_index >= 0 && position_index < static_cast<int>(position_names.size())) 
        {
            for (const auto& temp_set : temperatures) 
            {
                result.push_back(temp_set[position_index]);
            }
        }
        return result;
    }

    // 获取所有位置名称
    const std::vector<std::string>& getPositionNames() const { return position_names; }

    void setMeasurePoints(const std::vector<int>& positions, const std::vector<int>& indices) 
    {
        measure_positions_ = positions;
        measure_indices_ = indices;
    }

    const std::vector<int>& getMeasurePositions() const { return measure_positions_; }
    const std::vector<int>& getMeasureIndices() const { return measure_indices_; }
    int getMeasurePosition(int idx = 0) const { return idx < static_cast<int>(measure_positions_.size()) ? measure_positions_[idx] : 0; }
    int getMeasureIndex(int idx = 0) const { return idx < static_cast<int>(measure_indices_.size()) ? measure_indices_[idx] : 0; }
    const std::vector<double>& getTimes() const { return times; }
    size_t size() const { return times.size(); }
};

// 求解器接口
class heatTransferSolver
{
protected:

    // 时间步长
    double dt;
    int total_nodes;
    int total_steps;
    // 温度场
    std::vector<double> temperature;
    
    // 温度测量数据
    TemperatureMeasurement measurement;

    double initial_temp;

public:
    heatTransferSolver(double time_step)
    :
    dt(time_step),
    total_nodes(0),
    total_steps(0),
    initial_temp(25)
    {
        // 初始化完成，total_steps将在加载测量数据后设置
    }

    virtual ~heatTransferSolver() = default;

    virtual std::string typeName() const = 0;

    //传热过程求解
    virtual void solve() = 0;

    // 目标函数：计算模拟温度与测量温度的误差
    virtual void computeResiduals(std::vector<double>& residuals) = 0;

    virtual std::vector<double> rhoCp() const = 0;

    virtual std::vector<double> kappa() const = 0;
    
    // 设置初始温度
    void initialise(double iniT) 
    {
        initial_temp = iniT;
        temperature.resize(total_nodes, initial_temp);
    }

    bool loadMeasurementData
    (
        const std::string& filename, 
        const std::vector<int>& positions,
        const std::vector<int>& indices
    ) 
    {
        bool success = measurement.loadFromFile(filename);

        if (success) 
        {
            // 设置total_steps
            if (!measurement.getTimes().empty()) 
            {
                double total_time = measurement.getTimes().back();
                total_steps = std::round(total_time/dt);
            }
            
            measurement.setMeasurePoints(positions, indices);
        }
        return success;
    }
    
    inline void initialTemperature()
    {
        std::fill(temperature.begin(), temperature.end(), initial_temp);
    }

    // 获取温度结果
    inline const std::vector<double>& getTemperature() const 
    { 
        return temperature; 
    }
    
    inline double getInitialTemperature() const { return initial_temp; }

    inline double getTimeStep() const {return dt;}

    inline int getTotalNodes() const { return total_nodes; }

    inline int getTotalSteps() const { return total_steps; }

    inline const TemperatureMeasurement& getMeasurement() const { return measurement; }
};

// 自动微分兼容的求解器基类
template<typename T>
class heatTransferSolverAD 
{
protected:
    double dt;
    int total_nodes;
    int total_steps;
    int num_parameters;
    std::vector<T> temperature;
    const TemperatureMeasurement& measurement;
    T initial_temp;
    T zero;
    T one;

public:
    heatTransferSolverAD
    (
        heatTransferSolver& solver
    ) 
    : 
    dt(solver.getTimeStep()),
    total_nodes(solver.getTotalNodes()),
    total_steps(solver.getTotalSteps()),
    measurement(solver.getMeasurement()),
    initial_temp(solver.getInitialTemperature())
    {}

    virtual ~heatTransferSolverAD() = default;

    virtual void initialize(int num_parameters)
    {
        initial_temp.v.resize(num_parameters);
        initial_temp.v.setZero();
        zero.a = 0.0;
        zero.v.resize(num_parameters);
        zero.v.setZero();
        one.a = 1.0;
        one.v.resize(num_parameters);
        one.v.setZero();

        temperature.resize(total_nodes, initial_temp);
    }

    virtual void solve() = 0;
    
    // 自动微分兼容的残差
    virtual void computeResiduals(std::vector<T>& residuals) = 0;

    virtual void setupSolverParameters
    (
        const std::vector<T>& parameters,
        int param_index
    ) = 0;

    inline void initialTemperature()
    {
        std::fill(temperature.begin(), temperature.end(), initial_temp);
    }

    inline const TemperatureMeasurement& getMeasurement() const { return measurement; }
    
    // 获取时间步长
    inline T getTimeStep() const { return dt; }
    
    inline int getTotalNodes() const { return total_nodes; }

    inline int getTotalSteps() const { return total_steps; }
};

#endif // HEAT_TRANSFER_SOLVER_BASE_H
