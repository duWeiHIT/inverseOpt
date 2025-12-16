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

//常数
constexpr double pi(M_PI);
constexpr double piByTwo(M_PI_2);
constexpr double stefan_boltzmann(5.670367e-8);

// 边界条件类型枚举
enum class BoundaryType 
{
    FIXED_TEMPERATURE,           // 第一类：定壁温
    TIME_VARY_FIXED_TEMPERATURE, // 第二类：时变定壁温
    FIXED_HEAT_FLUX,             // 第三类：定热流
    CONVECTION,                  // 第四类：对流换热
    CONVECTION_RADIATION         // 第五类：对流换热+辐射
};

// 辐射模型枚举
enum class radiationModel
{
    NONE,
    DISCRETE_ORDINATE_METHOD,
    ROSSELAND
};

//前向声明
template<class T>
class radiationSolverBase;
class radiationSolver;
template<class T>
class radiationSolverAD;

// 边界条件结构体
struct BoundaryCondition 
{
    BoundaryType type;
    double value1;  // 温度或热流密度或流体温度
    double value2;  // 换热系数或辐射参数
    double value3;  // 辐射参数（发射率等)
};

// 材料属性类（分段线性）
class MaterialProperty 
{
private:
    std::vector<double> temperatures;
    std::vector<double> values; // 分段线性模式
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
    double refractiveIndex;                 // 折射率
    int    nodes;                           // 离散点数目
    MaterialProperty thermal_conductivity;  // 热导率
    MaterialProperty specific_heat;         // 比热
    MaterialProperty extinction;            // 衰减系数
    MaterialProperty albedo;                // 散射反照率系数
    std::string name;                       // 材料名称

    MaterialLayer() 
    :
    thickness(0.01),
    density(1000),
    refractiveIndex(1.0),
    nodes(10),
    thermal_conductivity("k_default"),
    specific_heat("cp_default"),
    extinction("extinction_default"),
    albedo("albedo_default"),
    name("default")
    {}
    
    MaterialLayer(const std::string& layer_name, double thick, int n, double rho) 
    :
    thickness(thick),
    density(rho),
    refractiveIndex(1.0),
    nodes(n),
    thermal_conductivity("k_" + layer_name),
    specific_heat("cp_" + layer_name),
    extinction("extinction_" + layer_name),
    albedo("albedo_" + layer_name),
    name(layer_name)
    {}

    void setRefractiveIndex(double n)
    {
        refractiveIndex = n;
    }
    
    inline double getK(double temp) const 
    {
        return thermal_conductivity.evaluate(temp);
    }
    
    inline double getCp(double temp) const 
    {
        return specific_heat.evaluate(temp);
    }

    inline double getExtinction(double temp) const 
    {
        return extinction.evaluate(temp);
    }

    inline double getAlbedo(double temp) const 
    {
        return albedo.evaluate(temp);
    }
};

// 自动微分兼容的材料层
template<typename T>
struct MaterialLayerAD 
{
    double thickness;
    double density;
    double refractiveIndex;
    int nodes;
    MaterialPropertyAD<T> thermal_conductivity;    // 热导率
    MaterialPropertyAD<T> specific_heat;          // 比热
    MaterialPropertyAD<T> extinction;      // 衰减系数
    MaterialPropertyAD<T> albedo;                // 散射反照率系数
    std::string name;

    MaterialLayerAD(const MaterialLayer& original_layer)
    :
    thickness(original_layer.thickness), 
    density(original_layer.density), 
    refractiveIndex(original_layer.refractiveIndex),
    nodes(original_layer.nodes),
    thermal_conductivity("k_" + original_layer.name),
    specific_heat("cp_" + original_layer.name),
    extinction("extinction_" + original_layer.name),
    albedo("albedo_" + original_layer.name),
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

        // 复制衰减系数温度点
        extinction.setData
        (
            original_layer.extinction.getTemperatures()
        );

        // 复制反照率温度点
        albedo.setData
        (
            original_layer.albedo.getTemperatures()
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

    inline T getExtinction(T temp) const 
    {
        return extinction.evaluate(temp);
    }

    inline T getAlbedo(T temp) const 
    {
        return albedo.evaluate(temp);
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

// 辐射强度类
template<typename T>
class intensity
{
    const radiationSolverBase<T>& solver_;
    double omega_; // 立体角
    double dAve_; // 平均辐射强度
    T zero_; // 零值
    T one_;  // 一值
    std::vector<T> I; // 辐射强度
    std::vector<T> A, B, C; //求解矩阵系数

public:
    intensity
    (
        const radiationSolverBase<T>& solver,
        const double dir,
        const double deltaPhi,
        const double deltaTheta
    )
    : solver_(solver),
    omega_(0.0),
    dAve_(0.0),
    zero_(T(0.0)),
    one_(T(1.0))
    {
        omega_ = 2.0*deltaPhi;
        // dAve_ = cosPhi*std::sin(0.5*deltaPhi)*(deltaTheta - std::cos(2.0*theta)*std::sin(deltaTheta));
        dAve_ = dir*deltaTheta;
    }

    void initialise(const T& zero, const T& one, int total_nodes)
    {
        zero_ = zero;
        one_ = one;
        I.resize(total_nodes, zero);
        A.resize(total_nodes, zero);
        B.resize(total_nodes, zero);
        C.resize(total_nodes, zero);
    }

    void applyBoundaryConditions(const std::vector<T>& temperature)
    {
        int n = I.size() - 1;
        if(dAve_ > 0)
        {
            //左边界
            const T& T0 = temperature[0];
            const double emissivity = solver_.getLeftEmissivity();
            I[0] = (emissivity*stefan_boltzmann/pi)*T0*T0*T0*T0
                    + ((1.0 - emissivity)/pi)*solver_.getIncidenceIntensityLeft();
            
            //右边界
            C[n] = zero_;
            A[n] = one_;
            B[n] = -one_;
        }
        else
        {
            //左边界
            C[0] = zero_;
            A[0] = one_;
            B[0] = -one_;

            //右边界
            const T& TN = temperature[n];
            const double emissivity = solver_.getRightEmissivity();
            I[n] = (emissivity*stefan_boltzmann/pi)*TN*TN*TN*TN
                    + ((1.0 - emissivity)/pi)*solver_.getIncidenceIntensityRight();
        }
    }

    void solve()
    {        
        // 计算辐射强度
        const std::vector<T>& temperature = solver_.getTemperature();
        
        for(size_t i = 1; i < I.size() - 1; ++i)
        {
            const T& Ti = temperature[i];
            const T beta = solver_.getExtinction(i);
            const T kappa = (one_ - solver_.getAlbedo(i))*beta;

            double dx = solver_.getDx(i);
            const double L = dAve_/dx;
            A[i] = L + beta*omega_;

            B[i] = -L*one_;
            C[i] = (1.0/pi)*omega_*stefan_boltzmann*kappa*Ti*Ti*Ti*Ti;
        }

        // 边界节点处理
        applyBoundaryConditions(temperature);

        if(dAve_ > 0)
        {
            for (size_t i = 1; i < I.size(); i++) 
            {
                I[i] = (C[i] - B[i]*I[i-1])/A[i];
            }
        }
        else
        {
            int n = I.size() - 1;
            for (int i = n-1; i >= 0; i--) 
            {
                I[i] = (C[i] - B[i]*I[i+1])/A[i];
            }
        }
    }

    const std::vector<T>& getIntensity() const { return I; }

    double omega() const { return omega_; }

    double dAve() const { return dAve_; }
};

// 辐射求解器基类
template<typename T>
class radiationSolverBase
{
protected:
    int total_nodes;
    const std::vector<T>& temperature_;
    double epsilon_left_;
    double epsilon_right_;
    std::vector<intensity<T>> Is_;

public:
    radiationSolverBase(const std::vector<T>& temperature)
    :total_nodes(temperature.size()),
    temperature_(temperature),
    epsilon_left_(1.0),
    epsilon_right_(1.0)
    {}

    inline const std::vector<T>& getTemperature() const
    {
        return temperature_;
    };

    inline const T& getIncidenceIntensityLeft() const
    {
        return Is_[1].getIntensity()[0];
    }

    inline const T& getIncidenceIntensityRight() const
    {
        return Is_[0].getIntensity()[total_nodes-1];
    };

    virtual double getLeftEmissivity() const = 0;

    virtual double getRightEmissivity() const = 0;

    virtual double getDx(int index) const = 0;

    virtual T getExtinction(int index) const = 0;

    virtual T getAlbedo(int index) const = 0;
};

// 求解器接口
class heatTransferSolver
{
protected:

    // 时间步长
    double dt;
    int total_nodes;
    int total_steps;
    radiationModel model;
    
    // 材料层
    std::vector<MaterialLayer> layers;

    // 温度场
    std::vector<double> temperature;
    
    // 温度测量数据
    TemperatureMeasurement measurement;

    std::unique_ptr<radiationSolver> radiation_solver;

    double initial_temp;

public:
    heatTransferSolver(double time_step, radiationModel model_type)
    :
    dt(time_step),
    total_nodes(0),
    total_steps(0),
    model(model_type),
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

    virtual double getLeftEmissivity() const = 0;

    virtual double getRightEmissivity() const = 0;

    virtual double getDx(int index = 0) const = 0;

    virtual double getExtinction(int index) const = 0;

    virtual double getAlbedo(int index) const = 0;
    
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
    
    inline radiationModel getRadiationModel() const
    {
        return model;
    }

    inline void initialTemperature()
    {
        std::fill(temperature.begin(), temperature.end(), initial_temp);
    }

    inline const std::vector<MaterialLayer>& getLayers() const
    {
        return layers;
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

// 辐射求解器
class radiationSolver
:public radiationSolverBase<double>
{
    const heatTransferSolver& solver_;

public:
    using radiationSolverBase<double>::Is_;

    radiationSolver(const heatTransferSolver& solver)
    : radiationSolverBase<double>(solver.getTemperature()),
        solver_(solver)
    {}

    void initialise()
    {
        const double deltaTheta = pi;
        const double deltaPhi = pi;
        Is_.push_back(intensity<double>(*this, 1.0, deltaPhi, deltaTheta));
        Is_.push_back(intensity<double>(*this, -1.0, deltaPhi, deltaTheta));
        double zero = 0.0;
        double one = 1.0;
        Is_[0].initialise(zero, one, this->total_nodes);
        Is_[1].initialise(zero, one, this->total_nodes);
    }

    virtual double getLeftEmissivity() const override
    {
        return solver_.getLeftEmissivity();
    }

    virtual double getRightEmissivity() const override
    {
        return solver_.getRightEmissivity();
    }

    virtual double getDx(int index) const override
    {
        return solver_.getDx(index);
    }
    
    virtual double getExtinction(int index) const override
    {
        return solver_.getExtinction(index);
    }

    virtual double getAlbedo(int index) const override
    {
        return solver_.getAlbedo(index);
    }

    void solve()
    {
        for(auto& I : Is_)
        {
            I.solve();
        }
    }

    void correct
    (
        std::vector<double>& D, 
        std::vector<double>& S,
        const std::vector<double>& temperature
    )
    {
        for(size_t i = 0; i < D.size(); ++i)
        {
            const double kappa = (1.0 - getAlbedo(i))*getExtinction(i);
            for(size_t j = 0; j < Is_.size(); ++j)
            {
                S[i] += Is_[j].getIntensity()[i]*kappa*Is_[j].omega();
            }
            D[i] += 4.0*stefan_boltzmann*kappa*std::pow(temperature[i],3);
        }
    }
};

// 自动微分兼容的求解器基类
template<typename T>
class heatTransferSolverAD 
{
protected:
    double dt;
    int total_nodes;
    int total_steps;
    radiationModel model;
    std::vector<MaterialLayerAD<T>> layers_ad;
    std::vector<T> temperature;
    const TemperatureMeasurement& measurement;
    std::unique_ptr<radiationSolverAD<T>> radiation_solver;
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
    model(solver.getRadiationModel()),
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

    virtual double getLeftEmissivity() const = 0;

    virtual double getRightEmissivity() const = 0;

    virtual double getDx(int index = 0) const = 0;

    virtual T getExtinction(int index) const = 0;

    virtual T getAlbedo(int index) const = 0;

    inline const std::vector<T>& getTemperature() const
    {
        return temperature;
    }

    inline std::vector<T>& getTemperature()
    {
        return temperature;
    }

    inline void initialTemperature()
    {
        std::fill(temperature.begin(), temperature.end(), initial_temp);
    }

    inline const std::vector<MaterialLayerAD<T>>& getLayers() const
    {
        return layers_ad;
    }

    inline const TemperatureMeasurement& getMeasurement() const { return measurement; }
    
    // 获取时间步长
    inline double getTimeStep() const { return dt; }
    
    inline int getTotalNodes() const { return total_nodes; }

    inline int getTotalSteps() const { return total_steps; }
};

//辐射求解器AD版本
template<typename T>
class radiationSolverAD
:public radiationSolverBase<T>
{
    const heatTransferSolverAD<T>& solver_;

public:
    using radiationSolverBase<T>::Is_;

    radiationSolverAD(const heatTransferSolverAD<T>& solver)
    :   radiationSolverBase<T>(solver.getTemperature()),
        solver_(solver)
    {}

    void initialise(int num_parameters)
    {
        const double deltaTheta = pi;
        const double deltaPhi = pi;
        
        Is_.push_back(intensity<T>(*this, 1.0, deltaPhi, deltaTheta));
        Is_.push_back(intensity<T>(*this, -1.0, deltaPhi, deltaTheta));

        T zero, one;
        zero.a = 0.0;
        zero.v.resize(num_parameters);
        zero.v.setZero();
        one.a = 1.0;
        one.v.resize(num_parameters);
        one.v.setZero();
        Is_[0].initialise(zero, one, this->total_nodes);
        Is_[1].initialise(zero, one, this->total_nodes);
    }

    virtual double getLeftEmissivity() const override
    {
        return solver_.getLeftEmissivity();
    }

    virtual double getRightEmissivity() const override
    {
        return solver_.getRightEmissivity();
    }

    virtual double getDx(int i) const override
    {
        return solver_.getDx(i);
    }
    
    virtual T getExtinction(int index) const override
    {
        return solver_.getExtinction(index);
    }

    virtual T getAlbedo(int index) const override
    {
        return solver_.getAlbedo(index);
    }

    void solve()
    {
        for(auto& I : Is_)
        {
            I.solve();
        }
    }

    void correct
    (
        std::vector<T>& D, 
        std::vector<T>& S,
        const std::vector<T>& temperature
    )
    {
        for(size_t i = 0; i < D.size(); ++i)
        {
            const T kappa = (1.0 - getAlbedo(i))*getExtinction(i);
            for(size_t j = 0; j < Is_.size(); ++j)
            {
                S[i] += Is_[j].getIntensity()[i]*kappa*Is_[j].omega();
            }
            D[i] += 4.0*stefan_boltzmann*kappa*ceres::pow(temperature[i],3);
        }
    }
};

#endif // HEAT_TRANSFER_SOLVER_BASE_H
