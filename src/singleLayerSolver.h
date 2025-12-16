#ifndef SINGLELAYER_SOLVER_H
#define SINGLELAYER_SOLVER_H

#include "heatTransferSolverBase.h"

// 一维非稳态导热求解器
class singleLayerSolver1D
:
public heatTransferSolver
{
private:
    // 几何和网格参数
    double length;
    
    //空间步长
    double dx;
    
    // 边界条件
    BoundaryCondition& left_boundary;
    BoundaryCondition& right_boundary;
    
    //求解矩阵系数
    std::vector<double> L, D, U, S;

    // 界面热导率计算（调和平均）
    inline double getInterfaceConductivity
    (
        double temp_left, 
        double temp_right
    ) const 
    {
        double k_left =  layers[0].getK(temp_left);
        double k_right = layers[0].getK(temp_right);
        if(model == radiationModel::ROSSELAND)
        {
            double n = layers[0].refractiveIndex;
            k_left +=
                16.0*n*n*stefan_boltzmann*pow(temp_left + 273.15,3)/
                (3.0*layers[0].getExtinction(temp_left));
            k_right += 
                16.0*n*n*stefan_boltzmann*pow(temp_right + 273.15,3)/
                (3.0*layers[0].getExtinction(temp_right));
        }
        return 2.0*k_left*k_right /(k_left + k_right);
    }
    
public:
    singleLayerSolver1D
    (
        const MaterialLayer& mat_layer,
        radiationModel radModel,
        double time_step,
        BoundaryCondition& left_bc, 
        BoundaryCondition& right_bc
    )
    : 
    heatTransferSolver(time_step, radModel),
    length(mat_layer.thickness),
    dx(0.0),
    left_boundary(left_bc), 
    right_boundary(right_bc)
    {
        layers.push_back(mat_layer);
        total_nodes = layers[0].nodes;
        dx = length/(total_nodes - 1);
    }

    virtual ~singleLayerSolver1D() = default;

    virtual std::string typeName() const override
    {
        return "singleLayerSolver1D";
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
        D.resize(total_nodes, 0.0);
        U.resize(total_nodes, 0.0);
        S.resize(total_nodes, 0.0);
    }
    
    // 边界条件处理函数
    inline void applyLeftBoundary
    (
        const std::vector<double>& T_old
    ) 
    {     
        const double rhoCp = layers[0].density*layers[0].getCp(T_old[0]);
        const double diag = 0.5*rhoCp/dt;   
        switch (left_boundary.type) 
        {
            case BoundaryType::FIXED_TEMPERATURE: 
            case BoundaryType::TIME_VARY_FIXED_TEMPERATURE:
            {
                D[0] = diag;
                U[0] = 0.0;
                S[0] = diag*left_boundary.value1;
                break;
            }
            case BoundaryType::FIXED_HEAT_FLUX: 
            {
                const double q = left_boundary.value1;
                const double k_0 = layers[0].getK(T_old[0]);
                const double k_1 = layers[0].getK(T_old[1]);
                const double k_avg = 2.0*k_0*k_1/(k_0 + k_1);

                const double upper = (k_avg+k_0)/(dx*dx);
                D[0] = diag + upper;
                U[0] = upper;
                S[0] = 2.0*q/dx + diag*T_old[0];
                break;
            }
            case BoundaryType::CONVECTION: 
            {
                const double h = left_boundary.value2;
                const double T_inf = left_boundary.value1;
                const double k_0 = layers[0].getK(T_old[0]);
                const double k_1 = layers[0].getK(T_old[1]);
                const double k_avg = 2.0*k_0*k_1/(k_0 + k_1);

                const double u1 = (k_avg+k_0)/(dx*dx);
                const double u2 = 2.0*h/dx;
                D[0] = diag + u1 + u2;
                U[0] = u1;
                S[0] = u2*T_inf + diag*T_old[0];
                break;
            }
            case BoundaryType::CONVECTION_RADIATION: 
            {
                const double h = left_boundary.value2;
                const double T_inf = left_boundary.value1;
                const double epsilon = left_boundary.value3;
                
                const double k_0 = layers[0].getK(T_old[0]);
                const double k_1 = layers[0].getK(T_old[1]);
                const double k_avg = 2.0*k_0*k_1/(k_0 + k_1);
                
                const double T_old_abs = T_old[0] + 273.15; 
                const double T_inf_abs = T_inf + 273.15;
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
        const double rhoCp = layers[0].density*layers[0].getCp(T_old[n]);
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
                const double k_n = layers[0].getK(T_old[n]);
                const double k_n_1 = layers[0].getK(T_old[n-1]);
                const double k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);

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

                const double k_n = layers[0].getK(T_old[n]);
                const double k_n_1 = layers[0].getK(T_old[n-1]);
                const double k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);

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
                
                const double k_n = layers[0].getK(T_old[n]);
                const double k_n_1 = layers[0].getK(T_old[n-1]);
                const double k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);
                
                const double T_old_abs = T_old[n] + 273.15; 
                const double T_inf_abs = T_inf + 273.15;

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
    
    inline void applyBoundaryCondition
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
        int n = total_nodes;
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
    
    // 每个时间步，主求解函数
    virtual void solve() override
    {           
        // 内部节点方程
        for (int i = 1; i < total_nodes - 1; i++) 
        {
            const double k_left = getInterfaceConductivity(temperature[i-1], temperature[i]);
            const double k_right = getInterfaceConductivity(temperature[i], temperature[i+1]);

            const double rhoCp = layers[0].density*layers[0].getCp(temperature[i]);
            const double diag = rhoCp/dt;
            const double r1 = 1.0/(dx*dx);
            const double u_left = k_left*r1;
            const double u_right = k_right*r1;
            
            L[i] = u_left;
            D[i] = diag + u_left + u_right;
            U[i] = u_right;
            S[i] = diag*temperature[i];
        }
            
        // 应用边界条件
        applyBoundaryCondition(temperature);
        
        // 添加辐射
        correctRad(D, S, temperature);
            
        // 求解方程组
        solveTDMA(temperature);
    }
    
    // 运行正演模型计算残差
    virtual void computeResiduals(std::vector<double>& residuals) override
    {
        initialTemperature();

        const auto& meas_times = measurement.getTimes();
        const auto& positions = measurement.getMeasurePositions();
    
        size_t meas_index = 0;
        size_t res_index = 0;
        bool time_vary_leftBC = (left_boundary.type == BoundaryType::TIME_VARY_FIXED_TEMPERATURE);
        bool time_vary_rightBC = (right_boundary.type == BoundaryType::TIME_VARY_FIXED_TEMPERATURE);
        
        for (int step = 0; step <= total_steps; ++step) 
        {
            double current_time = step*dt;
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
        return dx;
    }

    virtual double getExtinction(int index) const override
    {
        return layers[0].getExtinction(temperature[index]);
    }

    virtual double getAlbedo(int index) const override
    {
        return layers[0].getAlbedo(temperature[index]);
    }

    const MaterialLayer& getLayer() const { return layers[0]; }

    MaterialLayer& getLayer() { return layers[0]; }

    BoundaryCondition& getLeftBoundary() { return left_boundary; }

    BoundaryCondition& getRightBoundary() { return right_boundary; }

    double getLength() const { return length; }

    void printTemperatureProfile() const 
    {
        std::cout << "温度分布 (K):" << std::endl;
        for (int i = 0; i < total_nodes; i++) 
        {
            std::cout << "节点 " << i << " (x=" << i*dx << "m): " 
                      << temperature[i] << " K" << std::endl;
        }
    }
};

template<typename T>
class singleLayerSolver1DAD
:
public heatTransferSolverAD<T>
{
private:

    using heatTransferSolverAD<T>::dt;
    using heatTransferSolverAD<T>::total_nodes;
    using heatTransferSolverAD<T>::temperature;
    using heatTransferSolverAD<T>::layers_ad;
    // 几何和网格参数
    double length;
    
    //空间步长
    double dx;
    
    // 材料层
    const MaterialLayer& layer;
    
    // 边界条件
    BoundaryCondition& left_boundary;
    BoundaryCondition& right_boundary;
    

    //求解矩阵系数
    std::vector<T> L, D, U, S;

    // 界面热导率计算（调和平均）
    T getInterfaceConductivity
    (
        T temp_left, 
        T temp_right
    ) const 
    {
        T k_left = layers_ad[0].getK(temp_left);
        T k_right = layers_ad[0].getK(temp_right);
        if(this->model == radiationModel::ROSSELAND)
        {
            double n = layers_ad[0].refractiveIndex;
            k_left +=
                (16.0*n*n*stefan_boltzmann)*ceres::pow(temp_left + 273.15,3)/
                (3.0*layers_ad[0].getExtinction(temp_left));
            k_right += 
                (16.0*n*n*stefan_boltzmann)*ceres::pow(temp_right + 273.15,3)/
                (3.0*layers_ad[0].getExtinction(temp_right));
        }
        // 调和平均，确保热流连续性
        return 2.0*k_left*k_right/(k_left + k_right);
    }
    
public:
    singleLayerSolver1DAD
    (
        singleLayerSolver1D& original_solver
    )
    : 
    heatTransferSolverAD<T>(original_solver),
    length(original_solver.getLength()),
    dx(original_solver.getDx(0)),
    layer(original_solver.getLayer()),
    left_boundary(original_solver.getLeftBoundary()),
    right_boundary(original_solver.getRightBoundary())
    {}
    
    virtual ~singleLayerSolver1DAD() = default;

    virtual void initialize(int num_parameters)
    {
        heatTransferSolverAD<T>::initialize(num_parameters);
        
        if(this->model == radiationModel::DISCRETE_ORDINATE_METHOD)
        {
            this->radiation_solver = std::make_unique<radiationSolverAD<T>>(*this);
            this->radiation_solver->initialise(num_parameters);
        }

        // 初始化材料属性的AD变量
        const std::vector<double>& k_vals_original 
                = layer.thermal_conductivity.getValues();
        const std::vector<double>& cp_vals_original 
                = layer.specific_heat.getValues();

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

        layers_ad[0].thermal_conductivity.setValue(k_vals);
        layers_ad[0].specific_heat.setValue(cp_vals);

        if(!std::is_same<T, double>::value)
        {
            L.resize(total_nodes, this->zero);
            D.resize(total_nodes, this->zero);
            U.resize(total_nodes, this->zero);
            S.resize(total_nodes, this->zero);
        }
    }

    // 边界条件处理函数
    inline void applyLeftBoundary
    (
        const std::vector<T>& T_old
    ) 
    {        
        const T rhoCp = layers_ad[0].density*layers_ad[0].getK(T_old[0]);
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
                const T k_0 = layers_ad[0].getK(T_old[0]);
                const T k_1 = layers_ad[0].getK(T_old[1]);
                const T k_avg = 2.0*k_0*k_1/(k_0 + k_1);

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
                const T k_0 = layers_ad[0].getK(T_old[0]);
                const T k_1 = layers_ad[0].getK(T_old[1]);
                const T k_avg = 2.0*k_0*k_1/(k_0 + k_1);

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
                
                const T k_0 = layers_ad[0].getK(T_old[0]);
                const T k_1 = layers_ad[0].getK(T_old[1]);
                const T k_avg = 2.0*k_0*k_1/(k_0 + k_1);
                
                // 转换为绝对温度
                const T T_old_abs = T_old[0] + 273.15; 
                const double T_inf_abs = T_inf + 273.15;
                // 线性化辐射换热系数
                const T T_old_abs3 = T_old_abs*T_old_abs*T_old_abs;
                const T h_rad = 4.0*epsilon*stefan_boltzmann*T_old_abs3;
                const T source = 
                    epsilon*stefan_boltzmann*(pow(T_inf_abs, 4) - T_old_abs3*T_old_abs);

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
        const T rhoCp = layers_ad[0].density*layers_ad[0].getCp(T_old[n]);
        const T diag = 0.5*rhoCp/dt;    
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
                const T k_n = layers_ad[0].getK(T_old[n]);
                const T k_n_1 = layers_ad[0].getK(T_old[n-1]);
                const T k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);

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

                const T k_n = layers_ad[0].getK(T_old[n]);
                const T k_n_1 = layers_ad[0].getK(T_old[n-1]);
                const T k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);

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
                
                const T k_n = layers_ad[0].getK(T_old[n]);
                const T k_n_1 = layers_ad[0].getK(T_old[n-1]);
                const T k_avg = 2.0*k_n*k_n_1/(k_n + k_n_1);
                
                // 转换为绝对温度
                const T T_old_abs = T_old[n] + 273.15; 
                const double T_inf_abs = T_inf + 273.15;

                // 线性化辐射换热系数
                const T T_old_abs3 = T_old_abs*T_old_abs*T_old_abs;
                const T h_rad = 4.0*epsilon*stefan_boltzmann*T_old_abs3;
                const T source = 
                    epsilon*stefan_boltzmann*(pow(T_inf_abs, 4) - T_old_abs3*T_old_abs);

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
    
    inline void applyBoundaryCondition
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

    // TDMA算法求解三对角方程组
    inline void solveTDMA(std::vector<T>& x) 
    {
        int n = total_nodes;
        std::vector<T> P(n);
        std::vector<T> Q(n);
        
        // 前向消元
        P[0] = U[0]/D[0];
        Q[0] = S[0]/D[0];
        for (int i = 1; i < n; i++) 
        {
            const T denominator = D[i] - L[i]*P[i-1];
            
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
    
    // 每个时间步，主求解函数
    virtual void solve() override
    {         
        // 内部节点方程
        for (int i = 1; i < total_nodes - 1; i++) 
        {
            const T k_left = getInterfaceConductivity(temperature[i-1], temperature[i]);
            const T k_right = getInterfaceConductivity(temperature[i], temperature[i+1]);
            const T rhoCp = layers_ad[0].density*layers_ad[0].getCp(temperature[i]);
            const T diag = rhoCp/dt;
            const double rdx2 = 1.0/(dx*dx);
            const T u_left = k_left*rdx2;
            const T u_right = k_right*rdx2;
            
            L[i] = u_left;
            D[i] = diag + u_left + u_right;
            U[i] = u_right;
            S[i] = diag*temperature[i];
        }
        // 应用边界条件
        applyBoundaryCondition(temperature);
        
        // 添加辐射
        correctRad(D, S, temperature);

        // 求解方程组
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
        return dx;
    }

    virtual T getExtinction(int index) const override
    {
        return layers_ad[0].getExtinction(temperature[index]);
    }

    virtual T getAlbedo(int index) const override
    {
        return layers_ad[0].getAlbedo(temperature[index]);
    }

    // 运行正演模型计算残差
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
            // 保存上一个时间步的温度
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
            std::vector<T>& k_vals = layers_ad[0].thermal_conductivity.getValues();
        
            size_t k_count = k_vals.size();
            for (size_t i = 0; i < k_count; ++i) 
            {
                k_vals[i] = parameters[i];
            }
        }
        else if(param_index == 1)
        {
            std::vector<T>& beta_vals = layers_ad[0].extinction.getValues();
            for (size_t i = 0; i < beta_vals.size(); ++i) 
            {
                beta_vals[i] = parameters[i];
            }
        }
    }
};

#endif // SINGLELAYER_SOLVER_H
