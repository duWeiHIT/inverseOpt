
#ifndef MULTILAYER_SOLVER_H
#define MULTILAYER_SOLVER_H

#include "heatTransferSolverBase.h"

class multiLayerSolver1D
: public heatTransferSolver
{
private:

    double total_length;

    int sampleI;
  
    // 节点位置
    std::vector<double> node_positions; 
    // 节点所属层ID
    std::vector<int> node_layer_id; 
    
    // 边界条件
    BoundaryCondition& left_boundary;
    BoundaryCondition& right_boundary;
    
    // 求解矩阵系数
    std::vector<double> L, D, U, S;

    void initializeNodePositions() 
    {
        node_positions.resize(total_nodes);
        node_layer_id.resize(total_nodes);
        
        int offset = 0;
        int index = 0;
        node_positions[0] = 0.0;
        node_layer_id[0] = 0;
        for(const auto& layer:layers) 
        {
            const double dx = layer.thickness/layer.nodes;

            for(int i = 1; i < layer.nodes; i++)
            {
                node_positions[i + offset] += dx;
                node_layer_id[i + offset] = index++;
            }
            offset += layer.nodes - 1;
        }
    }
    
    // 获取节点处的材料属性
    inline double getThermalConductivity
    (
        int node_index, 
        double temp,
        int offset = 0
    ) const 
    {
        int layer_id = node_layer_id[node_index];
        double kappa = layers[layer_id + offset].getK(temp);
        
        if(model == radiationModel::ROSSELAND)
        {
            double n = layers[layer_id + offset].refractiveIndex;
            kappa +=
                16.0*n*n*stefan_boltzmann*pow(temp + 273.15,3)/
                (3.0*layers[layer_id + offset].getExtinction(temp));
        }

        return kappa;
    }
    
    inline double getSpecificHeat
    (
        int node_index, 
        double temp,
        int offset = 0
    ) const 
    {
        int layer_id = node_layer_id[node_index];
        return layers[layer_id + offset].getCp(temp);
    }
    
    inline double getDensity(int node_index, int offset = 0) const 
    {
        int layer_id = node_layer_id[node_index];
        return layers[layer_id + offset].density;
    }
    
    // 界面热导率计算（调和平均）
    inline double getInterfaceConductivity
    (
        int left_node, 
        int right_node, 
        double temp_left, 
        double temp_right,
        int offset = 0
    ) const 
    {
        double k_left = getThermalConductivity(left_node, temp_left, offset);
        double k_right = getThermalConductivity(right_node, temp_right);
        
        // 调和平均，确保热流连续性
        return 2.0*k_left*k_right/(k_left + k_right);
    }
 
    // 边界条件处理函数
    inline void applyLeftBoundary
    (
        const std::vector<double>& T_old
    ) 
    {   
        const double rhoCp = getDensity(0)*getSpecificHeat(0, T_old[0]);
        const double diag = 0.5*rhoCp/dt;     
        switch (left_boundary.type) 
        {
            case BoundaryType::FIXED_TEMPERATURE:
            case BoundaryType::TIME_VARY_FIXED_TEMPERATURE: 
            {
                // 第一类边界条件：定壁温

                D[0] = diag;
                U[0] = 0.0;
                S[0] = diag*left_boundary.value1;

                break;
            }
            case BoundaryType::FIXED_HEAT_FLUX: 
            {
                // 第二类边界条件：定热流
                const double q = left_boundary.value1;
                const double k_0 = getThermalConductivity(0, T_old[0]);
                const double k_1 = getThermalConductivity(1, T_old[1]);
                const double k_avg = 2.0*k_0*k_1/(k_0 + k_1);
                const double dx = node_positions[1] - node_positions[0];

                const double upper = (k_avg+k_0)/(dx*dx);
                D[0] = diag + upper;
                U[0] = upper;
                S[0] = 2.0*q/dx + diag*T_old[0];

                break;
            }
            case BoundaryType::CONVECTION: 
            {
                // 第三类边界条件：对流换热
                const double h = left_boundary.value2;
                const double T_inf = left_boundary.value1;
                const double k_0 = getThermalConductivity(0, T_old[0]);
                const double k_1 = getThermalConductivity(1, T_old[1]);
                const double k_avg = 2.0*k_0*k_1/(k_0 + k_1);
                const double dx = node_positions[1] - node_positions[0];

                const double u1 = (k_avg+k_0)/(dx*dx);
                const double u2 = 2.0*h/dx;
                D[0] = diag + u1 + u2;
                U[0] = u1;
                S[0] = u2*T_inf + diag*T_old[0];
                break;
            }
            case BoundaryType::CONVECTION_RADIATION: 
            {
                // 第四类边界条件：对流+辐射
                const double h = left_boundary.value2;
                const double T_inf = left_boundary.value1;
                const double epsilon = left_boundary.value3;
                
                const double k_0 = getThermalConductivity(0, T_old[0]);
                const double k_1 = getThermalConductivity(1, T_old[1]);
                const double k_avg = 2.0*k_0*k_1/(k_0 + k_1);
                const double dx = node_positions[1] - node_positions[0];
                
                // 转换为绝对温度
                const double T_old_abs = T_old[0] + 273.15; 
                const double T_inf_abs = T_inf + 273.15;
                // 线性化辐射换热系数
                const double T_old_abs3 = pow(T_old_abs, 3.0);
                const double h_rad = 4.0*epsilon*stefan_boltzmann*T_old_abs3;
                const double source = 
                    epsilon*stefan_boltzmann*(pow(T_inf_abs, 4) - T_old_abs3*T_old_abs);

                const double r1 = 2.0/dx;
                const double u1 = (k_avg+k_0)/(dx*dx);
                const double u2 = h*r1;
                const double u3 = h_rad*r1;
                
                D[0] = diag + u1 + u2 + u3;
                U[0] = u1;
                S[0] = u2*T_inf + u3*T_old_abs + r1*source + diag*T_old[0];
                break;
            }
        }
    }
    
    inline void applyRightBoundary
    (
        const std::vector<double>& T_old
    ) 
    {  
        int n = total_nodes - 1; 
        const double rhoCp = getDensity(n)*getSpecificHeat(n, T_old[n]);
        const double diag = 0.5*rhoCp/dt;       
        switch (right_boundary.type) 
        {
            case BoundaryType::FIXED_TEMPERATURE:
            case BoundaryType::TIME_VARY_FIXED_TEMPERATURE: 
            {
                L[n] = 0.0;
                D[n] = diag;
                S[n] = diag*right_boundary.value1;
                break;
            }
            case BoundaryType::FIXED_HEAT_FLUX: 
            {
                const double q = right_boundary.value1;
                const double k_n = getThermalConductivity(n, T_old[n]);
                const double k_n_1 = getThermalConductivity(n-1, T_old[n-1]);
                const double k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);
                const double dx = node_positions[n] - node_positions[n-1];

                const double upper = (k_avg+k_n_1)/(dx*dx);


                L[n] = upper;
                D[n] = diag + upper;
                S[n] = 2.0*q/dx + diag*T_old[n];
                break;
            }
            case BoundaryType::CONVECTION: 
            {
                const double h = right_boundary.value2;
                const double T_inf = right_boundary.value1;

                const double k_n = getThermalConductivity(n, T_old[n]);
                const double k_n_1 = getThermalConductivity(n-1, T_old[n-1]);
                const double k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);
                const double dx = node_positions[n] - node_positions[n-1];

                const double u1 = (k_avg+k_n)/(dx*dx);
                const double u2 = 2.0*h/dx; 
                
                L[n] = u1;
                D[n] = diag + u1 + u2;
                S[n] = u2*T_inf + diag*T_old[n];
                break;
            }
            case BoundaryType::CONVECTION_RADIATION: 
            {
                const double h = right_boundary.value2;
                const double T_inf = right_boundary.value1;
                const double epsilon = right_boundary.value3;
                
                const double k_n = getThermalConductivity(n, T_old[n]);
                const double k_n_1 = getThermalConductivity(n-1, T_old[n-1]);
                const double k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);
                const double dx = node_positions[n] - node_positions[n-1];
                
                // 转换为绝对温度
                const double T_old_abs = T_old[n] + 273.15; 
                const double T_inf_abs = T_inf + 273.15;

                // 线性化辐射换热系数
                const double T_old_abs3 = pow(T_old_abs, 3);
                const double h_rad = 4.0*epsilon*stefan_boltzmann*T_old_abs3;
                const double source = 
                    epsilon*stefan_boltzmann*(pow(T_inf_abs, 4.0) - T_old_abs3*T_old_abs);

                const double r1 = 2.0/dx;
                const double u1 = (k_avg+k_n_1)/(dx*dx);
                const double u2 = h*r1;
                const double u3 = h_rad*r1;

                L[n] = u1;
                D[n] = diag + u1 + u2 + u3;
                S[n] = u2*T_inf + u3*T_old_abs + r1*source + diag*T_old[n];
                break;
            }
        }
    }
    
    inline void applyBoundaryConditions
    (
        const std::vector<double>& T_old
    )
    {
        applyLeftBoundary(T_old);
        applyRightBoundary(T_old);
    }
    
    inline void correctRad
    (
        std::vector<double>& D, 
        std::vector<double>& S,
        const std::vector<double>& temperature
    )
    {
        if(model == radiationModel::DISCRETE_ORDINATE_METHOD)
        {
            radiation_solver->solve();
            radiation_solver->correct(D, S, temperature);
        }
    }

    // TDMA算法求解三对角方程组
    inline void solveTDMA(std::vector<double>& x) 
    {
        const int n = total_nodes;
        std::vector<double> P(n);
        std::vector<double> Q(n);
        
        // 前向消元
        P[0] = U[0]/D[0];
        Q[0] = S[0]/D[0];
        
        for (int i = 1; i < n; i++) 
        {
            const double denominator = D[i] - L[i]*P[i-1];
            P[i] = U[i]/denominator;
            Q[i] = (S[i] + L[i]*Q[i-1])/denominator;
        }
        
        // 回代
        x[n-1] = Q[n-1];
        for (int i = n-2; i >= 0; i--) 
        {
            x[i] = P[i]*x[i+1] + Q[i];
        }
    }
  
public:

    multiLayerSolver1D
    (
        const std::vector<MaterialLayer>& material_layers,
        radiationModel modle_type,
        double time_step,
        BoundaryCondition& left_bc,
        BoundaryCondition& right_bc
    ) 
    : 
    heatTransferSolver(time_step, modle_type),
    total_length(0.0),
    sampleI(0),
    left_boundary(left_bc), 
    right_boundary(right_bc) 
    {
        //赋值构造
        layers = material_layers;

        // 计算总长度和总节点数目
        total_length = 0.0;
        for(const auto& layer : layers) 
        {
            total_length += layer.thickness;
            total_nodes += layer.nodes;
        }
        // 移除界面处重复的节点
        total_nodes -= layers.size() - 1;
                
        // 初始化节点位置和层ID
        initializeNodePositions();
    }

    virtual ~multiLayerSolver1D() = default;

    virtual std::string typeName() const override
    {
        return "multiLayerSolver1D";
    }

    //初始化
    void initialise(double iniT) 
    {
        heatTransferSolver::initialise(iniT);
        
        if(model == radiationModel::DISCRETE_ORDINATE_METHOD)
        {
            radiation_solver = std::make_unique<radiationSolver>(*this);
            radiation_solver->initialise();
        }

        // 初始化矩阵和温度场
        L.resize(total_nodes, 0.0);
        U.resize(total_nodes, 0.0);
        D.resize(total_nodes, 0.0);
        S.resize(total_nodes, 0.0);
    }
  
    // 主求解函数（考虑耦合界面处理）
    virtual void solve() override
    {
        temperature = temperature;
        
        for (int i = 1; i < total_nodes - 1; i++) 
        {
            const double dx_left = node_positions[i] - node_positions[i-1];
            const double dx_right = node_positions[i+1] - node_positions[i];
            const double dx_avg = 0.5*(dx_left + dx_right);

            // 使用调和平均
            const double k_eff_left = getInterfaceConductivity(i-1, i, temperature[i-1], temperature[i]);
            double k_eff_right;

            double rhoCp = getDensity(i)*getSpecificHeat(i, temperature[i]);
            
            // 检查是否为界面节点
            bool is_interface = (node_layer_id[i+1] != node_layer_id[i]);
            
            //节点位于界面，控制容积内的比热和密度由两侧介质计算
            if (is_interface) 
            {
                const double s1 = 0.5*dx_left/dx_avg;
                const double s2 = 0.5*dx_right/dx_avg;
                rhoCp *= s1;
                const double rhoRight = getDensity(i, 1);
                const double cpRight = getSpecificHeat(i, temperature[i], 1);
                rhoCp += rhoRight*cpRight*s2;
                k_eff_right = getInterfaceConductivity(i, i+1, temperature[i], temperature[i+1], 1);
            }
            else
            {
                k_eff_right = getInterfaceConductivity(i, i+1, temperature[i], temperature[i+1]);
            }

            const double diag = rhoCp/dt;
            const double u_left = k_eff_left/(dx_left*dx_avg);
            const double u_right = k_eff_right/(dx_right*dx_avg);
            
            L[i] = u_left;
            D[i] = diag + u_left + u_right;
            U[i] = u_right;
            S[i] = diag*temperature[i];
        }
        
        // 应用边界条件
        applyBoundaryConditions(temperature);

        //添加辐射
        correctRad(D, S, temperature);
        
        // 求解
        solveTDMA(temperature);
    }
    
    // 运行正演模型计算残差
    virtual void computeResiduals(std::vector<double>& residuals) override
    {
        initialTemperature();

        const auto& meas_times = this->measurement.getTimes();
        const auto& positions = this->measurement.getMeasurePositions();
        
        size_t meas_index = 0;
        size_t res_index = 0;

        bool time_vary_leftBC = (left_boundary.type == BoundaryType::TIME_VARY_FIXED_TEMPERATURE);
        bool time_vary_rightBC = (right_boundary.type == BoundaryType::TIME_VARY_FIXED_TEMPERATURE);

        for (int step = 0; step <= total_steps; ++step) 
        {
            double current_time = step * dt;
            
            if(time_vary_leftBC)
            {
                left_boundary.value1 = measurement.getTemperatureByTime(current_time, 0);
            }

            if(time_vary_rightBC)
            {
                int last_pos = measurement.getPositionNames().size() - 1;
                right_boundary.value1 = measurement.getTemperatureByTime(current_time, last_pos);
            }

            solve();
            
            // 在测量时间点计算残差
            if(meas_index < meas_times.size() && 
                std::abs(current_time - meas_times[meas_index]) < dt/2.0) 
            {
                for (size_t p = 0; p < positions.size(); ++p) 
                {
                    int pos_idx = measurement.getMeasureIndices()[p];
                    residuals[res_index++] = 
                        temperature[positions[p]] - measurement.getTemperature(meas_index, pos_idx);
                }
                meas_index++;
            }
        }
    }
    
    virtual double getLeftEmissivity() const override
    {
        return left_boundary.value3;
    }

    virtual double getRightEmissivity() const override
    {
        return right_boundary.value3;
    }

    virtual double getDx(int index) const override
    {
        return 0.5*(node_positions[index + 1] - node_positions[index - 1]);
    }

    virtual double getExtinction(int index) const override
    {
        int layer_id = node_layer_id[index];
        int layer_id1 = node_layer_id[index + 1];
        if(layer_id != layer_id1)
        {
            // 界面节点，取两侧介质的加权平均
            const double dx_left = node_positions[index] - node_positions[index - 1];
            const double dx_right = node_positions[index + 1] - node_positions[index];
            const double dx_avg = 0.5*(dx_left + dx_right);
            const double s1 = 0.5*dx_left/dx_avg;
            const double s2 = 0.5*dx_right/dx_avg;
            double ext_left = layers[layer_id].getExtinction(temperature[index]);
            double ext_right = layers[layer_id1].getExtinction(temperature[index]);
            return ext_left*s1 + ext_right*s2;
        }
        return layers[layer_id].getExtinction(temperature[index]);
    }

    virtual double getAlbedo(int index) const override
    {
        int layer_id = node_layer_id[index];
        int layer_id1 = node_layer_id[index + 1];
        if(layer_id != layer_id1)
        {
            // 界面节点，取两侧介质的加权平均
            const double dx_left = node_positions[index] - node_positions[index - 1];
            const double dx_right = node_positions[index + 1] - node_positions[index];
            const double dx_avg = 0.5*(dx_left + dx_right);
            const double s1 = 0.5*dx_left/dx_avg;
            const double s2 = 0.5*dx_right/dx_avg;
            double albedo_left = layers[layer_id].getAlbedo(temperature[index]);
            double albedo_right = layers[layer_id1].getAlbedo(temperature[index]);
            return albedo_left*s1 + albedo_right*s2;
        }

        return layers[layer_id].getAlbedo(temperature[index]);
    }

    // 获取材料层信息
    const MaterialLayer& getLayer() const { return layers[sampleI]; }

    MaterialLayer& getLayer() { return layers[sampleI]; }

    BoundaryCondition& getLeftBoundary() { return left_boundary; }
    BoundaryCondition& getRightBoundary() { return right_boundary; }
    const BoundaryCondition& getLeftBoundary() const{ return left_boundary; }
    const BoundaryCondition& getRightBoundary() const{ return right_boundary; }
    const std::vector<MaterialLayer>& getLayers() const { return layers; }

    void setSampleIndex(int i){sampleI = i;}

    int getSampleIndex(){return sampleI;}

    const std::vector<double>& getNodePositions() const {return node_positions;}

    const std::vector<int>& getNodeLayerId() const {return node_layer_id;}

    double getTotalLength() const { return total_length; }
    
    // 打印多层温度分布
    void printTemperatureProfile() const 
    {
        std::cout << "多层结构温度分布 (K):" << std::endl;
        for (int i = 0; i < total_nodes; i++) 
        {
            std::cout << "节点 " << i << " (x=" << node_positions[i] << "m, 层" 
                      << node_layer_id[i] << "): " << temperature[i] << " K" << std::endl;
        }
    }
    
    // 打印层信息
    void printLayerInfo() const 
    {
        std::cout << "多层结构信息:" << std::endl;
        for (size_t i = 0; i < layers.size(); ++i) 
        {
            std::cout << "层 " << i << " [" << layers[i].name << "]: " 
                      << "厚度=" << layers[i].thickness << "m, "
                      << "节点数=" << layers[i].nodes << std::endl;
        }
        std::cout << "总长度: " << total_length << "m, 总节点数: " << total_nodes << std::endl;
    }
};

// 自动微分兼容的求解器
template<typename T>
class multiLayerSolver1DAD 
: public heatTransferSolverAD<T>
{
private:
    using heatTransferSolverAD<T>::dt;
    using heatTransferSolverAD<T>::total_nodes;
    using heatTransferSolverAD<T>::temperature;
    using heatTransferSolverAD<T>::layers_ad;

    double total_length;
    int sampleI = 0;
    const std::vector<MaterialLayer>& layers;
    const std::vector<double>& node_positions; 
    const std::vector<int>& node_layer_id; 
    
    BoundaryCondition& left_boundary;
    BoundaryCondition& right_boundary;
    
    std::vector<T> L, D, U, S;

    // 自动微分兼容的材料属性获取
    T getThermalConductivity(int node_index, T temp, int offset = 0) const 
    {
        int layer_id = node_layer_id[node_index] + offset;
        T kappa = layers_ad[layer_id].thermal_conductivity.evaluate(temp);
        if(this->model == radiationModel::ROSSELAND)
        {
            double n = layers_ad[layer_id + offset].refractiveIndex;
            kappa +=
                16.0*n*n*stefan_boltzmann*ceres::pow(temp + 273.15,3)/
                (3.0*layers_ad[layer_id + offset].getExtinction(temp));
        }
        return kappa;
    }
    
    T getSpecificHeat(int node_index, T temp, int offset = 0) const 
    {
        int layer_id = node_layer_id[node_index] + offset;
        return layers_ad[layer_id].specific_heat.evaluate(temp);
    }
    
    double getDensity(int node_index, int offset = 0) const 
    {
        int layer_id = node_layer_id[node_index] + offset;
        return layers_ad[layer_id].density;
    }
    
    T getInterfaceConductivity
    (
        int left_node, 
        int right_node, 
        T temp_left,
        T temp_right,
        int offset = 0
    ) const 
    {
        T k_left = getThermalConductivity(left_node, temp_left, offset);
        T k_right = getThermalConductivity(right_node, temp_right);
        
        return 2.0*k_left*k_right/(k_left + k_right);
    }

    // 边界条件处理函数
    inline void applyLeftBoundary
    (
        const std::vector<T>& T_old
    ) 
    {        
        const T rhoCp = getDensity(0)*getSpecificHeat(0, T_old[0]);
        const T diag = (0.5/dt)*rhoCp;
        switch (left_boundary.type) 
        {
            case BoundaryType::FIXED_TEMPERATURE:
            case BoundaryType::TIME_VARY_FIXED_TEMPERATURE: 
            {
                // 第一类边界条件：定壁温
                D[0] = diag;
                U[0] = this->zero;
                S[0] = diag*left_boundary.value1;

                break;
            }
            case BoundaryType::FIXED_HEAT_FLUX: 
            {
                // 第二类边界条件：定热流
                const double q = left_boundary.value1;
                const T k_0 = getThermalConductivity(0, T_old[0]);
                const T k_1 = getThermalConductivity(1, T_old[1]);
                const T k_avg = 2.0*k_0*k_1/(k_0 + k_1);
                const double dx = node_positions[1] - node_positions[0];
                const double rdx = 1.0/dx;

                const T upper = (k_avg+k_0)*(rdx*rdx);
                D[0] = diag + upper;
                U[0] = upper;
                S[0] = 2.0*q*rdx + diag*T_old[0];

                break;
            }
            case BoundaryType::CONVECTION: 
            {
                // 第三类边界条件：对流换热
                const double h = left_boundary.value2;
                const double T_inf = left_boundary.value1;
                const T k_0 = getThermalConductivity(0, T_old[0]);
                const T k_1 = getThermalConductivity(1, T_old[1]);
                const T k_avg = 2.0*k_0*k_1/(k_0 + k_1);
                const double dx = node_positions[1] - node_positions[0];
                const double rdx = 1.0/dx;

                const T u1 = (k_avg+k_0)*(rdx*rdx);
                const double u2 = 2.0*h*rdx;
                D[0] = diag + u1 + u2;
                U[0] = u1;
                S[0] = u2*T_inf + diag*T_old[0];
                break;
            }
            case BoundaryType::CONVECTION_RADIATION: 
            {
                // 第四类边界条件：对流+辐射
                const double h = left_boundary.value2;
                const double T_inf = left_boundary.value1;
                const double epsilon = left_boundary.value3;
                
                const T k_0 = getThermalConductivity(0, T_old[0]);
                const T k_1 = getThermalConductivity(1, T_old[1]);
                const T k_avg = 2.0*k_0*k_1/(k_0 + k_1);
                const double dx = node_positions[1] - node_positions[0];
                
                // 转换为绝对温度
                const T T_old_abs = T_old[0] + 273.15; 
                const double T_inf_abs = T_inf + 273.15;
                // 线性化辐射换热系数
                const T T_old_abs3 = T_old_abs*T_old_abs*T_old_abs;
                const T h_rad = (4.0*epsilon*stefan_boltzmann)*T_old_abs3;
                const T source = 
                    (epsilon*stefan_boltzmann)*(pow(T_inf_abs, 4) - T_old_abs3*T_old_abs);

                const double rdx = 2.0/dx;
                const T u1 = (k_avg+k_0)*(rdx*rdx*0.25);
                const double u2 = h*rdx;
                const T u3 = h_rad*rdx;
                
                D[0] = diag + u1 + u2 + u3;
                U[0] = u1;
                S[0] = u2*T_inf + u3*T_old_abs + rdx*source + diag*T_old[0];
                break;
            }
        }
    }
    
    inline void applyRightBoundary
    (
        const std::vector<T>& T_old
    ) 
    {  
        int n = total_nodes - 1;
        const T rhoCp = getDensity(n)*getSpecificHeat(n, T_old[n]);  
        const T diag = (0.5/dt)*rhoCp;      
        switch (right_boundary.type) 
        {
            case BoundaryType::FIXED_TEMPERATURE:
            case BoundaryType::TIME_VARY_FIXED_TEMPERATURE: 
            {
                L[n] = this->zero;
                D[n] = diag;
                S[n] = diag*right_boundary.value1;
                break;
            }
            case BoundaryType::FIXED_HEAT_FLUX: 
            {
                const double q = right_boundary.value1;
                const T k_n = getThermalConductivity(n, T_old[n]);
                const T k_n_1 = getThermalConductivity(n-1, T_old[n-1]);
                const T k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);
                const double dx = node_positions[n] - node_positions[n-1];
                const double rdx = 1.0/dx;

                const T upper = (k_avg+k_n_1)*(rdx*rdx);

                L[n] = upper;
                D[n] = diag + upper;
                S[n] = 2.0*q*rdx + diag*T_old[n];
                break;
            }
            case BoundaryType::CONVECTION: 
            {
                const double h = right_boundary.value2;
                const double T_inf = right_boundary.value1;

                const T k_n = getThermalConductivity(n, T_old[n]);
                const T k_n_1 = getThermalConductivity(n-1, T_old[n-1]);
                const T k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);
                const double dx = node_positions[n] - node_positions[n-1];
                const double rdx = 1.0/dx;

                const T u1 = (k_avg+k_n)*(rdx*rdx);
                const double u2 = 2.0*h*rdx; 
                
                L[n] = u1;
                D[n] = diag + u1 + u2;
                S[n] = u2*T_inf + diag*T_old[n];
                break;
            }
            case BoundaryType::CONVECTION_RADIATION: 
            {
                const double h = right_boundary.value2;
                const double T_inf = right_boundary.value1;
                const double epsilon = right_boundary.value3;
                
                const T k_n = getThermalConductivity(n, T_old[n]);
                const T k_n_1 = getThermalConductivity(n-1, T_old[n-1]);
                const T k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);

                const double dx = node_positions[n] - node_positions[n-1];
                
                // 转换为绝对温度
                const T T_old_abs = T_old[n] + 273.15; 
                const double T_inf_abs = T_inf + 273.15;

                // 线性化辐射换热系数
                const T T_old_abs3 = T_old_abs*T_old_abs*T_old_abs;
                const T h_rad = (4.0*epsilon*stefan_boltzmann)*T_old_abs3;
                const T source = 
                    (epsilon*stefan_boltzmann)*(pow(T_inf_abs, 4) - T_old_abs3*T_old_abs);

                const double rdx = 2.0/dx;
                const T u1 = (k_avg+k_n_1)*(rdx*rdx*0.25);
                const double u2 = h*rdx;
                const T u3 = h_rad*rdx;

                L[n] = u1;
                D[n] = diag + u1 + u2 + u3;
                S[n] = u2*T_inf + u3*T_old_abs + rdx*source + diag*T_old[n];
                break;
            }
        }
    }
    
    inline void applyBoundaryConditions
    (
        const std::vector<T>& T_old
    )
    {
        applyLeftBoundary(T_old);
        applyRightBoundary(T_old);
    }

    inline void correctRad
    (
        std::vector<T>& D, 
        std::vector<T>& S,
        const std::vector<T>& temperature
    )
    {
        if(this->model == radiationModel::DISCRETE_ORDINATE_METHOD)
        {
            this->radiation_solver->solve();
            this->radiation_solver->correct(D, S, temperature);
        }
    }

    // TDMA求解器（自动微分兼容）
    void solveTDMA(std::vector<T>& x) 
    {
        int n = total_nodes;
        std::vector<T> P(n), Q(n);
        
        P[0] = U[0]/D[0];
        Q[0] = S[0]/D[0];
        
        for (int i = 1; i < n; i++) 
        {
            const T denominator = D[i] - L[i]*P[i-1];
            P[i] = U[i] / denominator;
            Q[i] = (S[i] + L[i] * Q[i-1])/denominator;
        }
        
        x[n-1] = Q[n-1];
        for (int i = n-2; i >= 0; i--) 
        {
            x[i] = P[i]*x[i+1] + Q[i];
        }
    }

public:
    // 从原始求解器构造自动微分版本
    multiLayerSolver1DAD
    (
        multiLayerSolver1D& original_solver
    )
    : 
    heatTransferSolverAD<T>(original_solver),
    total_length(original_solver.getTotalLength()),
    sampleI(original_solver.getSampleIndex()),
    layers(original_solver.getLayers()),
    node_positions(original_solver.getNodePositions()),
    node_layer_id(original_solver.getNodeLayerId()),
    left_boundary(original_solver.getLeftBoundary()),
    right_boundary(original_solver.getRightBoundary())
    {}

    virtual ~multiLayerSolver1DAD() = default;

    virtual void initialize(int num_parameters)
    {
        heatTransferSolverAD<T>::initialize(num_parameters);
        
        if(this->model == radiationModel::DISCRETE_ORDINATE_METHOD)
        {
            this->radiation_solver = std::make_unique<radiationSolverAD<T>>(*this);
            this->radiation_solver->initialise(num_parameters);
        }
        // 初始化材料属性的AD变量
        for(size_t i = 0; i < layers.size(); ++i)
        {
            layers_ad.push_back(MaterialLayerAD<T>(layers[i]));
        }

        for(size_t i = 0; i < layers_ad.size(); ++i)
        {
            auto& layer = layers_ad[i];
            auto& original_layer = layers[i];
            const std::vector<double>& k_vals_original 
                = original_layer.thermal_conductivity.getValues();
            const std::vector<double>& cp_vals_original 
                = original_layer.specific_heat.getValues();

            std::vector<T> k_vals(k_vals_original.size());
            std::vector<T> cp_vals(cp_vals_original.size());
            if(std::is_same<T, double>::value)
            {
                for(size_t j = 0; j < k_vals_original.size(); ++j)
                {
                    k_vals[j] = T(k_vals_original[j]);
                }
                for(size_t j = 0; j < cp_vals_original.size(); ++j)
                {
                    cp_vals[j] = T(cp_vals_original[j]);
                }
            }
            else
            {   
                for(size_t j = 0; j < k_vals_original.size(); ++j)
                {
                    k_vals[j] = T(k_vals_original[j]);
                    k_vals[j].v.resize(num_parameters);
                    k_vals[j].v.setZero();
                }
                for(size_t j = 0; j < cp_vals_original.size(); ++j)
                {
                    cp_vals[j] = T(cp_vals_original[j]);
                    cp_vals[j].v.resize(num_parameters);
                    cp_vals[j].v.setZero();
                }
            }

            layer.thermal_conductivity.setValue(k_vals);
            layer.specific_heat.setValue(cp_vals);
        }

        if(!std::is_same<T, double>::value)
        {
            L.resize(total_nodes, this->zero);
            D.resize(total_nodes, this->zero);
            U.resize(total_nodes, this->zero);
            S.resize(total_nodes, this->zero);
        }
    }

    virtual void solve() override
    {        
        for (int i = 1; i < total_nodes - 1; i++) 
        {
            const double dx_left = node_positions[i] - node_positions[i-1];
            const double dx_right = node_positions[i+1] - node_positions[i];
            const double dx_avg = 0.5*(dx_left + dx_right);

            const T k_eff_left = getInterfaceConductivity(i-1, i, temperature[i-1], temperature[i]);
            T k_eff_right;

            T rhoCp = getDensity(i)*getSpecificHeat(i, temperature[i]);

            bool is_interface = (node_layer_id[i+1] != node_layer_id[i]);
            //节点位于界面，控制容积内的比热和密度由两侧介质计算
            if (is_interface) 
            {
                const double s1 = 0.5*dx_left/dx_avg;
                const double s2 = 0.5*dx_right/dx_avg;
                rhoCp *= s1;
                const double rhoRight = getDensity(i, 1);
                const T cpRight = getSpecificHeat(i, temperature[i], 1);
                rhoCp += (rhoRight*s2)*cpRight;
                k_eff_right = getInterfaceConductivity(i, i+1, temperature[i], temperature[i+1], 1);
            }
            else
            {
                k_eff_right = getInterfaceConductivity(i, i+1, temperature[i], temperature[i+1]);
            }
            
            const T diag = rhoCp/dt;
            const T u1 = k_eff_left/(dx_left*dx_avg);
            const T u2 = k_eff_right/(dx_right*dx_avg);
            
            L[i] = u1;
            D[i] = diag + u1 + u2;
            U[i] = u2;
            S[i] = diag*temperature[i];
        }
        
        // 应用边界条件
        applyBoundaryConditions(temperature);
        
        //添加辐射
        correctRad(D, S, temperature);
        
        solveTDMA(temperature);
    }

    virtual double getLeftEmissivity() const override
    {
        return left_boundary.value3;
    }

    virtual double getRightEmissivity() const override
    {
        return right_boundary.value3;
    }

    virtual double getDx(int index) const override
    {
        return 0.5*(node_positions[index + 1] - node_positions[index - 1]);
    }

    virtual T getExtinction(int index) const override
    {
        int layer_id = node_layer_id[index];
        int layer_id1 = node_layer_id[index + 1];
        if(layer_id != layer_id1)
        {
            // 界面节点，取两侧介质的加权平均
            const double dx_left = node_positions[index] - node_positions[index - 1];
            const double dx_right = node_positions[index + 1] - node_positions[index];
            const double dx_avg = 0.5*(dx_left + dx_right);
            const double s1 = 0.5*dx_left/dx_avg;
            const double s2 = 0.5*dx_right/dx_avg;
            T ext_left = layers_ad[layer_id].getExtinction(temperature[index]);
            T ext_right = layers_ad[layer_id1].getExtinction(temperature[index]);
            return ext_left*s1 + ext_right*s2;
        }
        return layers_ad[layer_id].getExtinction(temperature[index]);
    }

    virtual T getAlbedo(int index) const override
    {
        int layer_id = node_layer_id[index];
        int layer_id1 = node_layer_id[index + 1];
        if(layer_id != layer_id1)
        {
            // 界面节点，取两侧介质的加权平均
            const double dx_left = node_positions[index] - node_positions[index - 1];
            const double dx_right = node_positions[index + 1] - node_positions[index];
            const double dx_avg = 0.5*(dx_left + dx_right);
            const double s1 = 0.5*dx_left/dx_avg;
            const double s2 = 0.5*dx_right/dx_avg;
            T albedo_left = layers_ad[layer_id].getAlbedo(temperature[index]);
            T albedo_right = layers_ad[layer_id1].getAlbedo(temperature[index]);
            return albedo_left*s1 + albedo_right*s2;
        }

        return layers_ad[layer_id].getAlbedo(temperature[index]);
    }

    virtual void computeResiduals(std::vector<T>& residuals) override
    {
        this->initialTemperature();
        
        const auto& meas_times = this->measurement.getTimes();
        const auto& positions = this->measurement.getMeasurePositions();

        size_t meas_index = 0;
        size_t res_index = 0;

        bool time_vary_leftBC = (left_boundary.type == BoundaryType::TIME_VARY_FIXED_TEMPERATURE);
        bool time_vary_rightBC = (right_boundary.type == BoundaryType::TIME_VARY_FIXED_TEMPERATURE);
        
        for (int step = 0; step <= this->total_steps; ++step) 
        {
            double current_time = step*dt;
            if(time_vary_leftBC)
            {
                left_boundary.value1 = this->measurement.getTemperatureByTime(current_time, 0);
            }

            if(time_vary_rightBC)
            {
                int last_pos = this->measurement.getPositionNames().size() - 1;
                right_boundary.value1 = this->measurement.getTemperatureByTime(current_time, last_pos);
            }

            solve();
            
            // 在测量时间点计算残差
            if(meas_index < meas_times.size() && 
                std::abs(current_time - meas_times[meas_index]) < dt/2.0) 
            {
                for (size_t p = 0; p < positions.size(); ++p) 
                {
                    int pos_idx = this->measurement.getMeasureIndices()[p];
                    residuals[res_index++] = 
                        temperature[positions[p]] - this->measurement.getTemperature(meas_index, pos_idx);
                }
                meas_index++;
            }
        }
    }
    
    virtual void setupSolverParameters
    (
        const std::vector<T>& parameters,
        int param_index
    ) override
    {
        if(param_index == 0)
        {
            MaterialPropertyAD<T>& k_prop = layers_ad[sampleI].thermal_conductivity;
            std::vector<T>& k_vals = k_prop.getValues();
            
            size_t k_count = k_vals.size();
            for (size_t i = 0; i < k_count; ++i) 
            {
                k_vals[i] = parameters[i];
            }
        }
        else if(param_index == 1)
        {
            MaterialPropertyAD<T>& cp_prop = layers_ad[sampleI].specific_heat;
            std::vector<T>& cp_vals = cp_prop.getValues();
            for (size_t i = 0; i < cp_vals.size(); ++i) 
            {
                cp_vals[i] = parameters[i];
            }
        }
    }
};


#endif // MULTILAYER_SOLVER_H
