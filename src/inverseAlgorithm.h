#ifndef INVERSSE_ALGORITHM_H
#define INVERSSE_ALGORITHM_H

#include <memory>
#include <functional>
#include <unordered_map>
#include <random>
#include <Eigen/Dense>
#include "singleLayerSolver.h"
#include "multiLayerSolver.h"

// 热物性反演专用目标函数
class ThermalPropertyObjective
{
private:
    heatTransferSolver& solver_;
    MaterialLayer& layer;
    int residualSize;
    int inversion_mode_;
    double k_lower_bound_;
    double k_upper_bound_;
    double beta_lower_bound_;
    double beta_upper_bound_;
    double regularization_lambda_;
public:
    ThermalPropertyObjective
    (
        heatTransferSolver& solver,
        MaterialLayer& material_layer,
        int inversion_mode
    )
    : 
    solver_(solver),
    layer(material_layer),
    residualSize(solver.getMeasurement().getTimes().size() * solver.getMeasurement().getMeasurePositions().size()),
    inversion_mode_(inversion_mode),
    k_lower_bound_(1.0),
    k_upper_bound_(1000.0),
    beta_lower_bound_(100.0),
    beta_upper_bound_(50000.0),
    regularization_lambda_(0.001)
    {}
    
    // 设置参数边界
    void setThermalConductivityBounds(double lower, double upper) {
        k_lower_bound_ = lower;
        k_upper_bound_ = upper;
    }

    void setExtinctionBounds(double lower, double upper) {
        beta_lower_bound_ = lower;
        beta_upper_bound_ = upper;
    }
    
    void setRegularizationParameter(double lambda) {
        regularization_lambda_ = lambda;
    }
    
    // 获取参数边界
    double getThermalConductivityLowerBound() const { return k_lower_bound_; }
    double getThermalConductivityUpperBound() const { return k_upper_bound_; }
    double getExtinctionLowerBound() const { return beta_lower_bound_; }
    double getExtinctionUpperBound() const { return beta_upper_bound_; }
    
    // 应用参数边界约束
    std::vector<double> applyBounds(const std::vector<double>& parameters) const 
    {
        std::vector<double> bounded_params = parameters;
        MaterialProperty& k_prop_ = layer.thermal_conductivity;
        size_t k_count = k_prop_.getValues().size();
        
        // 约束热导率
        for (size_t i = 0; i < k_count && i < parameters.size(); ++i) 
        {
            bounded_params[i] = std::clamp(parameters[i], 0.0, 1.0);
        }
        
        // 约束消光系数
        for (size_t i = k_count; i < bounded_params.size(); ++i) 
        {
            bounded_params[i] = std::clamp(parameters[i], 0.0, 1.0);
        }

        return bounded_params;
    }

    // 计算残差向量
    std::vector<double> computeResiduals(const std::vector<double>& parameters) 
    {
        if (static_cast<int>(parameters.size()) != getParameterCount()) 
        {
            throw std::invalid_argument("Parameter count mismatch");
        }
        updateMaterialProperties(parameters);
        std::vector<double> residuals(residualSize);
        solver_.computeResiduals(residuals);
        return residuals;
    }

    void computeResiduals(const double* const* parameters, double* residuals) 
    {
        const MaterialProperty& k_prop_ = layer.thermal_conductivity;
        size_t k_count = k_prop_.getValues().size();
        std::vector<double> params_vec(getParameterCount());
        for (size_t i = 0; i < k_count; ++i) 
        {
            params_vec[i] = parameters[0][i];
        }

        if(inversion_mode_ == 0)
        {
            const MaterialProperty& beta_prop = layer.extinction;
            for (size_t i = 0; i < beta_prop.getValues().size(); ++i) 
            {
                params_vec[i + k_count] = parameters[1][i];
            }
        }

        const auto res = computeResiduals(params_vec);
        if(residuals)
        {
            for (size_t i = 0; i < res.size(); ++i) 
            {
                residuals[i] = res[i];
            }
        }
    }

    void updateMaterialProperties(const std::vector<double>& parameters) 
    {
        //应用边界约束
        const std::vector<double> bounded_params = applyBounds(parameters);
        MaterialProperty& k_prop = layer.thermal_conductivity;
        std::vector<double>& k_vals = k_prop.getValues();
        size_t k_count = k_vals.size();
        
        std::vector<double> vals(bounded_params.size());
        const double deltaK = k_upper_bound_ - k_lower_bound_;
        
        for(size_t i = 0; i < k_vals.size(); ++i)
        {
            vals[i] = bounded_params[i]*deltaK + k_lower_bound_;
        }

        const double deltaBeta = beta_upper_bound_ - beta_lower_bound_;
        for(size_t i = k_count; i < parameters.size(); ++i)
        {
            vals[i] = bounded_params[i]*deltaBeta + beta_lower_bound_;
        }

        for (size_t i = 0; i < k_count; ++i) 
        {
            k_vals[i] = vals[i];
        }

        std::vector<double>& beta_vals = layer.extinction.getValues(); 

        for (size_t i = 0; i < beta_vals.size() && i + k_count < bounded_params.size(); ++i) 
        {
            beta_vals[i] = vals[i+k_count];
        }
    }

    double evaluate(const std::vector<double>& parameters) 
    {
        std::vector<double> residuals = computeResiduals(parameters);
        
        double sum_squared = 0.0;
        for (const auto& r : residuals) 
        {
            sum_squared += r * r;
        }
        
        double regularization = 0.0;

        // 添加Tikhonov正则化（光滑性约束）
        if (regularization_lambda_ > 0) 
        {
            const MaterialProperty& k_prop = layer.thermal_conductivity;
            // 对热导率添加二阶差分正则化（鼓励光滑变化）
            const auto& k_vals = k_prop.getValues();
            for (size_t i = 1; i < k_vals.size() - 1; ++i) 
            {
                double second_diff = k_vals[i-1] - 2*k_vals[i] + k_vals[i+1];
                regularization += regularization_lambda_ * second_diff * second_diff;
            }

            if(inversion_mode_ == 0)
            {
                const MaterialProperty& beta_prop = layer.extinction;
                // 对消光系数添加二阶差分正则化（鼓励光滑变化）
                const auto& beta_vals = beta_prop.getValues();
                for (size_t i = 1; i < beta_vals.size() - 1; ++i) 
                {
                    double second_diff = beta_vals[i-1] - 2*beta_vals[i] + beta_vals[i+1];
                    regularization += regularization_lambda_ * second_diff * second_diff;
                }
            }
        }
        
        return sqrt(sum_squared / residuals.size()) + regularization;
    }

    int getParameterCount() const 
    {
        const MaterialProperty& k_prop = layer.thermal_conductivity;
        const MaterialProperty& beta_prop = layer.extinction;
        if (inversion_mode_ == 0) 
        {
           return k_prop.getValues().size() + beta_prop.getValues().size();
        }

        return  k_prop.getValues().size();;
    }

    heatTransferSolver& solver() { return solver_; }

    const MaterialProperty& kappa() const { return layer.thermal_conductivity; }

    const MaterialProperty& cp() const { return layer.specific_heat; }

    const MaterialProperty& extinction() const { return layer.extinction; }

    const MaterialProperty& albedo() const { return layer.albedo; }

    inline int getInversionMode() const { return inversion_mode_; }

    inline double getRegularizationParameter() const { return regularization_lambda_; }
    
    //归一化
    void normalizeParameters(std::vector<double>& parameters) const
    {
        const MaterialProperty& k_prop = layer.thermal_conductivity;

        const double deltaK = k_upper_bound_ - k_lower_bound_;
        for(size_t i = 0; i < k_prop.getValues().size(); ++i)
        {
            parameters[i] = (parameters[i] - k_lower_bound_)/deltaK;
        }

        if(inversion_mode_ == 0)
        {   
            size_t k_count = k_prop.getValues().size();
            const double deltaBeta = beta_upper_bound_ - beta_lower_bound_;
            for(size_t i = k_count; i < parameters.size(); ++i)
            {
                parameters[i] = (parameters[i] - beta_lower_bound_)/deltaBeta;
            }
        }
    }

    int getResidualCount() const
    {
        return residualSize;
    }
};

// 反演结果结构
struct InversionResult 
{
    std::vector<double> optimal_parameters;
    double optimal_value;
    int iterations;
    bool converged;
    std::string message;
    std::vector<double> history; // 收敛历史
};

// 算法配置
struct AlgorithmConfig 
{
    int max_iterations = 500;
    double tolerance = 1e-4;
    int population_size = 50;
    int maxReStart = 2;
    double mutation_rate = 0.7;
    double crossover_rate = 0.001;

    std::string algorithm_type = "LEVENBERG_MARQUARDT"; // 可选值见下
    std::string linear_solver = "DENSE_QR";            // 线性求解器选择
};

// 反演算法接口
class InversionAlgorithm 
{
protected:

    AlgorithmConfig config;

public:
    virtual ~InversionAlgorithm() = default;
    virtual InversionResult optimize
    (
        ThermalPropertyObjective* objective, 
        const std::vector<double>& initial_guess
    ) = 0;
    virtual std::string getName() const = 0;

    //输出算法配置
    virtual void outputConfig() const = 0;
    //配置算法
    virtual void setConfig() = 0;
};

// 遗传算法
class GeneticAlgorithm : public InversionAlgorithm 
{
private:
    std::default_random_engine generator_;
    std::uniform_real_distribution<double> distribution_;
    
    // 个体结构
    struct Individual 
    {
        std::vector<double> parameters;
        double fitness;
        
        bool operator<(const Individual& other) const 
        {
            return fitness < other.fitness;
        }
    };
    
    // 初始化种群
    std::vector<Individual> initializePopulation
    (
        int population_size, 
        int param_count,
        ThermalPropertyObjective* objective
    ) 
    {
        std::vector<Individual> population(population_size);
        
        // 获取边界信息
        for (int i = 0; i < population_size; ++i) 
        {
            population[i].parameters.resize(param_count);
              
            size_t k_count = objective->kappa().getTemperatures().size();
                
            // 初始化热导率参数
            std::uniform_real_distribution<double> k_dist(0.0, 1.0);
            for (size_t j = 0; j < k_count && j < static_cast<size_t>(param_count); ++j) 
            {
                population[i].parameters[j] = k_dist(generator_);
            }
                
            // 初始化比热参数
            std::uniform_real_distribution<double> cp_dist(0.0, 1.0);
            for (size_t j = k_count; j < static_cast<size_t>(param_count); ++j) 
            {
                population[i].parameters[j] = cp_dist(generator_);
            }

            //从小到大排序
            std::sort(population[i].parameters.begin(), population[i].parameters.end());
            
            population[i].fitness = objective->evaluate(population[i].parameters);
        }
        
        return population;
    }
    
    // 选择操作（锦标赛选择）
    Individual tournamentSelection(const std::vector<Individual>& population, int tournament_size) 
    {
        std::uniform_int_distribution<int> index_dist(0, population.size() - 1);
        Individual best = population[index_dist(generator_)];
        
        for (int i = 1; i < tournament_size; ++i) 
        {
            Individual candidate = population[index_dist(generator_)];
            if (candidate.fitness < best.fitness) 
            { 
                // 最小化问题
                best = candidate;
            }
        }
        return best;
    }
    
    // 交叉操作（模拟二进制交叉）
    std::vector<double> simulatedBinaryCrossover
    (
        const std::vector<double>& parent1,
        const std::vector<double>& parent2,
        double crossover_rate, double eta = 2.0
    ) 
    {
        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
        std::uniform_real_distribution<double> u_dist(0.0, 1.0);
        
        std::vector<double> child(parent1.size());
        
        for (size_t i = 0; i < parent1.size(); ++i) 
        {
            if (prob_dist(generator_) < crossover_rate) 
            {
                double u = u_dist(generator_);
                double beta;
                if (u <= 0.5) 
                {
                    beta = std::pow(2.0 * u, 1.0 / (eta + 1.0));
                } else 
                {
                    beta = std::pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta + 1.0));
                }
                
                child[i] = 0.5 * ((1 + beta) * parent1[i] + (1 - beta) * parent2[i]);
            } 
            else 
            {
                child[i] = parent1[i];
            }
        }
        return child;
    }
    
    // 多项式变异
    void polynomialMutation
    (
        std::vector<double>& individual, 
        double mutation_rate, 
        double mutation_scale = 0.1, 
        double eta_m = 20.0
    ) 
    {
        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
        std::uniform_real_distribution<double> u_dist(0.0, 1.0);
        
        for (size_t i = 0; i < individual.size(); ++i) 
        {
            if (prob_dist(generator_) < mutation_rate) 
            {
                double u = u_dist(generator_);
                double delta;
                
                if (u < 0.5) 
                {
                    delta = std::pow(2.0 * u, 1.0 / (eta_m + 1.0)) - 1.0;
                } 
                else 
                {
                    delta = 1.0 - std::pow(2.0 * (1.0 - u), 1.0 / (eta_m + 1.0));
                }
                
                individual[i] += mutation_scale * delta;
            }
        }
    }

public:
    GeneticAlgorithm() : distribution_(0.0, 1.0) 
    {
        generator_.seed(std::random_device{}());
    }
    
    InversionResult optimize
    (
        ThermalPropertyObjective* objective,
        const std::vector<double>& initial_guess
    ) override 
    {
        
        InversionResult result;
        int n = initial_guess.size();
        int population_size = config.population_size;
        
        // 初始化种群
        auto population = initializePopulation(population_size, n, objective);
        
        // 记录最佳个体
        Individual global_best = population[0];
        for (const auto& ind : population) 
        {
            if (ind.fitness < global_best.fitness) 
            {
                global_best = ind;
            }
        }
        
        result.history.push_back(global_best.fitness);
        
        for (int iter = 0; iter < config.max_iterations; ++iter) 
        {
            std::vector<Individual> new_population;
            
            // 精英保留：直接保留最佳个体
            new_population.push_back(global_best);
            
            // 生成新一代
            while (static_cast<int>(new_population.size()) < population_size) 
            {
                // 选择
                Individual parent1 = tournamentSelection(population, 3);
                Individual parent2 = tournamentSelection(population, 3);
                    
                // 交叉
                auto child_params = 
                    simulatedBinaryCrossover(parent1.parameters, parent2.parameters, config.crossover_rate);
                    
                // 变异
                polynomialMutation(child_params, config.mutation_rate);
                
                // 应用边界约束
                child_params = objective->applyBounds(child_params);
                    
                Individual child;
                child.parameters = child_params;
                child.fitness = objective->evaluate(child.parameters);
                    
                new_population.push_back(child);
            }
            
            // 更新种群
            population = new_population;
            
            // 更新全局最佳
            for (const auto& ind : population) 
            {
                if (ind.fitness < global_best.fitness) 
                {
                    global_best = ind;
                }
            }
            
            result.history.push_back(global_best.fitness);
            
            // 检查收敛
            if (iter > 10) 
            {
                double improvement = result.history[iter-10] - global_best.fitness;
                if (std::abs(improvement) < config.tolerance) 
                {
                    result.optimal_parameters = global_best.parameters;
                    result.optimal_value = global_best.fitness;
                    result.iterations = iter + 1;
                    result.converged = true;
                    result.message = "收敛于迭代 " + std::to_string(iter + 1);
                    return result;
                }
            }
        }
        
        result.optimal_parameters = global_best.parameters;
        result.optimal_value = global_best.fitness;
        result.iterations = config.max_iterations;
        result.converged = false;
        result.message = "达到最大迭代次数";
        return result;
    }
    
    virtual void outputConfig() const override
    {
        std::cout << "遗传算法配置：" << std::endl;
        std::cout << "  收敛残差： " << config.tolerance << std::endl;
        std::cout << "  最大迭代次数: " << config.max_iterations << std::endl;
        std::cout << "  种群大小: " << config.population_size << std::endl;
        std::cout << "  变异率: " << config.mutation_rate << std::endl;
        std::cout << "  交叉率: " << config.crossover_rate << std::endl;
    }

    virtual void setConfig() override
    {
        std::cout << "输入遗传算法参数"<<std::endl;
        std::cout << "输入种群大小: "<< std::endl;
        std::cin >> config.population_size;
        std::cout << "输入最大迭代次数: "<< std::endl;
        std::cin >> config.max_iterations;
        std::cout << "输入变异率: "<< std::endl;
        std::cin >> config.mutation_rate;
        std::cout << "输入交叉率: "<< std::endl;
        std::cin >> config.crossover_rate;
    }

    std::string getName() const override { return "遗传算法"; }
};

// 粒子群优化算法
class ParticleSwarmOptimization 
: public InversionAlgorithm 
{
private:
    std::default_random_engine generator_;
    
    // PSO参数
    double w; // 惯性权重
    double c1; // 个体学习因子
    double c2; // 社会学习因子

    struct Particle 
    {
        std::vector<double> position;
        std::vector<double> velocity;
        std::vector<double> best_position;
        double best_fitness;
        double current_fitness;
    };

public:
    
    ParticleSwarmOptimization()
    {
        generator_.seed(std::random_device{}());

        w = 0.729; // 惯性权重
        c1 = 1.494; // 个体学习因子
        c2 = 1.494; // 社会学习因子
    }

    InversionResult optimize
    (
        ThermalPropertyObjective* objective,
        const std::vector<double>& initial_guess
    ) override
    {
        InversionResult result;
        int n = initial_guess.size();
        int population_size = config.population_size;
        
        std::vector<Particle> particles(population_size);
        std::vector<double> global_best_position = initial_guess;
        double global_best_fitness = objective->evaluate(initial_guess);
        
        // 初始化粒子
        std::uniform_real_distribution<double> pos_dist(-10.0, 10.0);
        std::uniform_real_distribution<double> vel_dist(-1.0, 1.0);
        
        for (int i = 0; i < population_size; ++i) 
        {
            particles[i].position.resize(n);
            particles[i].velocity.resize(n);
            
            for (int j = 0; j < n; ++j) 
            {
                particles[i].position[j] = pos_dist(generator_);
                particles[i].velocity[j] = vel_dist(generator_);
            }
            
            particles[i].best_position = particles[i].position;
            particles[i].current_fitness = objective->evaluate(particles[i].position);
            particles[i].best_fitness = particles[i].current_fitness;
            
            if (particles[i].best_fitness < global_best_fitness) 
            {
                global_best_fitness = particles[i].best_fitness;
                global_best_position = particles[i].position;
            }
        }
        
        result.history.push_back(global_best_fitness);
        
        std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
        
        for (int iter = 0; iter < config.max_iterations; ++iter) 
        {
            for (auto& particle : particles) 
            {
                // 更新速度
                for (int j = 0; j < n; ++j) 
                {
                    double r1 = prob_dist(generator_);
                    double r2 = prob_dist(generator_);
                    
                    particle.velocity[j] = w * particle.velocity[j] +
                                          c1 * r1 * (particle.best_position[j] - particle.position[j]) +
                                          c2 * r2 * (global_best_position[j] - particle.position[j]);
                }
                
                // 更新位置
                for (int j = 0; j < n; ++j) 
                {
                    particle.position[j] += particle.velocity[j];
                }
                
                // 评估新位置
                particle.current_fitness = objective->evaluate(particle.position);
                
                // 更新个体最佳
                if (particle.current_fitness < particle.best_fitness) 
                {
                    particle.best_fitness = particle.current_fitness;
                    particle.best_position = particle.position;
                    
                    // 更新全局最佳
                    if (particle.best_fitness < global_best_fitness) {
                        global_best_fitness = particle.best_fitness;
                        global_best_position = particle.position;
                    }
                }
            }
            
            result.history.push_back(global_best_fitness);
            
            // 检查收敛
            if (iter > 10) 
            {
                double improvement = result.history[iter-10] - global_best_fitness;
                if (std::abs(improvement) < config.tolerance) 
                {
                    result.optimal_parameters = global_best_position;
                    result.optimal_value = global_best_fitness;
                    result.iterations = iter + 1;
                    result.converged = true;
                    result.message = "收敛于迭代 " + std::to_string(iter + 1);
                    return result;
                }
            }
        }
        
        result.optimal_parameters = global_best_position;
        result.optimal_value = global_best_fitness;
        result.iterations = config.max_iterations;
        result.converged = false;
        result.message = "达到最大迭代次数";
        return result;
    }

    virtual void outputConfig() const override
    {
        std::cout << "粒子群优化算法配置：" << std::endl;
        std::cout << "  收敛残差： " << config.tolerance << std::endl;
        std::cout << "  最大迭代次数: " << config.max_iterations << std::endl;
        std::cout << "  种群大小: " << config.population_size << std::endl;
        std::cout << "  惯性权重: " << w << std::endl;
        std::cout << "  个体学习因子: " << c1 << std::endl;
        std::cout << "  社会学习因子: " << c2 << std::endl;
    }

    virtual void setConfig() override
    {
        std::cout << "输入粒子群优化算法参数"<<std::endl;
        std::cout << "输入种群大小: "<< std::endl;
        std::cin >> config.population_size;
        std::cout << "输入最大迭代次数: "<< std::endl;
        std::cin >> config.max_iterations;
        std::cout << "输入惯性权重: "<< std::endl;
        std::cin >> w;
        std::cout << "输入个体学习因子: "<< std::endl;
        std::cin >> c1;
        std::cout << "输入社会学习因子: "<< std::endl;
        std::cin >> c2;
    }
    
    std::string getName() const override { return "粒子群优化算法"; }
};

//CMA-ES算法
class CMAES : public InversionAlgorithm 
{
private:
    std::default_random_engine generator_;
    int population_size_;
    double sigma_;  // 初始步长
    
    // CMA-ES核心参数
    double mu_eff_;     // 有效选择强度
    double c_sigma_;    // 步长学习率
    double d_sigma_;    // 步长阻尼参数
    double cc_;         // 进化路径学习率
    double c1_;         // 秩1更新学习率
    double cmu_;        // 秩μ更新学习率
    
    Eigen::VectorXd mean_;          // 分布均值
    Eigen::MatrixXd covariance_;    // 协方差矩阵
    Eigen::MatrixXd B_;              // 特征向量矩阵
    Eigen::VectorXd D_;              // 特征值平方根对角线
    Eigen::VectorXd pc_;             // 进化路径（协方差）
    Eigen::VectorXd p_sigma_;        // 进化路径（步长）
    
    // 重启机制相关参数
    int max_restarts_;               // 最大重启次数
    int current_restart_;           // 当前重启次数
    double initial_sigma_;          // 初始步长（保存原始值）
    int initial_population_size_;    // 初始种群大小
    Eigen::VectorXd best_mean_;     // 历史最佳均值
    double best_fitness_;           // 历史最佳适应度
    int no_improvement_count_;      // 无改善代数计数
    int evaluations_since_restart_; // 当前重启后的评估次数
    
    // 特征分解相关
    bool eigen_decomp_updated_;
    
public:
    CMAES(int population_size = 0, double sigma = 0.1) 
        : population_size_(population_size), sigma_(sigma), max_restarts_(1),
          current_restart_(0), initial_sigma_(sigma), initial_population_size_(population_size),
          best_fitness_(std::numeric_limits<double>::max()), no_improvement_count_(0),
          evaluations_since_restart_(0), eigen_decomp_updated_(false)
    {
        max_restarts_ = config.maxReStart;
        generator_.seed(std::random_device{}());
    }
    
    InversionResult optimize(ThermalPropertyObjective* objective, const std::vector<double>& initial_guess) override 
    {
        InversionResult result;
        int n = initial_guess.size();
        
        // 自动设置种群大小（CMA-ES推荐值）
        if (population_size_ == 0) 
        {
            population_size_ = 4 + std::floor(3 * std::log(n));
            initial_population_size_ = population_size_;
        }
        
        // 初始化最佳解记录
        best_mean_ = Eigen::VectorXd::Map(initial_guess.data(), n);
        best_fitness_ = std::numeric_limits<double>::max();
        current_restart_ = 0;
        no_improvement_count_ = 0;
        
        std::cout << "CMA-ES开始优化..." << std::endl;
        std::cout << "参数维度: " << n << ", 最大重启次数: " << max_restarts_ << std::endl;
        
        // 重启循环
        while (current_restart_ <= max_restarts_) 
        {
            if(current_restart_ > 0)
            {
                std::cout << "=== 第 " << current_restart_ << " 次重启 ===" << std::endl;

                // 设置当前重启的参数
                setupRestartParameters();
            }
            
            // 初始化CMA-ES参数
            initializeParameters(n, current_restart_ == 0 ? 
                std::vector<double>(best_mean_.data(), best_mean_.data() + n) : 
                perturbInitialGuess(objective, best_mean_, n));
            
            // 设置权重
            Eigen::VectorXd weights = setupWeights(population_size_);
            mu_eff_ = 1.0 / (weights.array().square().sum());
            
            // 基于mu_eff调整参数
            adjustParametersBasedOnDimension(n);
            
            evaluations_since_restart_ = 0;
            
            // 代数循环
            bool converged = false;
            for (int gen = 0; gen < config.max_iterations && !converged; ++gen) 
            {
                // 生成新种群
                std::vector<Eigen::VectorXd> population;
                std::vector<double> fitness;
                
                generatePopulation(population, fitness, objective, n);
                evaluations_since_restart_ += population_size_;
                
                // 选择最佳个体更新分布
                updateDistribution(population, fitness, weights);
                
                // 记录最佳适应度
                double current_best_fitness = *std::min_element(fitness.begin(), fitness.end());
                result.history.push_back(current_best_fitness);
                
                // 更新历史最佳解
                if (current_best_fitness < best_fitness_) 
                {
                    best_fitness_ = current_best_fitness;
                    best_mean_ = mean_;
                    no_improvement_count_ = 0;
                } 
                else 
                {
                    no_improvement_count_++;
                }
                
                // 检查收敛条件
                converged = checkConvergence(gen, current_best_fitness);
                
                if (gen % 20 == 0 || converged) 
                {
                    std::cout << "代 " << gen << ": 当前最佳 = " << current_best_fitness 
                             << ", 历史最佳 = " << best_fitness_ 
                             << ", 步长 = " << sigma_ << std::endl;
                }
                
                // 检查是否需要提前重启（陷入局部最优）
                if (shouldEarlyRestart(gen, current_best_fitness)) 
                {
                    std::cout << "检测到局部最优，提前重启..." << std::endl;
                    break;
                }
            }
            
            // 检查总体收敛条件
            if (checkGlobalConvergence() == 1) 
            {
                std::cout << "全局收敛条件满足，优化结束。" << std::endl;
                break;
            }
            else if(checkGlobalConvergence() == -1)
            {
                std::cout << "未达到收敛残差，优化结束。" << std::endl;
                break;
            }
            else
            {
                if(current_restart_ <= max_restarts_)
                {
                    std::cout << "全局收敛条件未满足，继续优化..." << std::endl;
                }
                else
                {
                    std::cout << "达到最大重启次数，优化结束。" << std::endl;
                }
            }
            
            // 准备下一次重启
            current_restart_++;
            if (current_restart_ <= max_restarts_) 
            {
                std::cout << "准备第 " << current_restart_ << " 次重启..." << std::endl;
            }
        }
        
        // 获取最终结果（使用历史最佳解）
        std::vector<double> best_params(best_mean_.data(), best_mean_.data() + best_mean_.size());
        result.optimal_parameters = best_params;
        result.optimal_value = objective->evaluate(best_params);
        result.iterations = result.history.size();
        result.converged = checkGlobalConvergence() == 1;
        result.message = "CMA-ES算法完成（重启次数: " + std::to_string(current_restart_) + ")";
        
        std::cout << "优化完成，总重启次数: " << current_restart_ 
                 << ", 最终最佳适应度: " << best_fitness_ << std::endl;
        
        return result;
    }
    
private:
    void setupRestartParameters() 
    {
        // IPOP-CMA-ES策略：每次重启增加种群大小
        population_size_ = initial_population_size_ * std::pow(2, current_restart_);
        population_size_ = std::min(population_size_, 1000); // 防止过大
            
        // 重置步长，但可能稍微调整
        sigma_ = initial_sigma_ * std::pow(2.0, current_restart_ * 0.3);
        sigma_ = std::max(sigma_, 1e-3);

        std::cout << "重启 " << current_restart_ << ": 种群大小 = " << population_size_ 
                << ", 步长 = " << sigma_ << std::endl;
    }
    
    std::vector<double> perturbInitialGuess
    (
        ThermalPropertyObjective* objective,
        const Eigen::VectorXd& best_mean, 
        int n
    ) 
    {
        // 对最佳解加入随机扰动，帮助探索新区域
        std::vector<double> perturbed_guess(n);
        std::normal_distribution<double> normal(0.0, 0.1*sigma_);
        
        for (int i = 0; i < n; ++i) 
        {
            perturbed_guess[i] = best_mean(i) + normal(generator_);
        }

        return objective->applyBounds(perturbed_guess);
    }
    
    void initializeParameters(int n, const std::vector<double>& initial_guess) 
    {
        mean_ = Eigen::VectorXd::Map(initial_guess.data(), n);
        covariance_ = Eigen::MatrixXd::Identity(n, n);
        B_ = Eigen::MatrixXd::Identity(n, n);
        D_ = Eigen::VectorXd::Ones(n);
        pc_ = Eigen::VectorXd::Zero(n);
        p_sigma_ = Eigen::VectorXd::Zero(n);
        eigen_decomp_updated_ = true;
        
        // 重启时重置进化路径
        if (current_restart_ > 0) 
        {
            pc_.setZero();
            p_sigma_.setZero();
        }
    }
    
    void adjustParametersBasedOnDimension(int n) 
    {
        // 标准的CMA-ES参数设置
        c_sigma_ = (mu_eff_ + 2.0) / (n + mu_eff_ + 3.0);
        d_sigma_ = 1.0 + c_sigma_ + 2.0*std::max(0.0, std::sqrt((mu_eff_ - 1.0) / (n + 1.0)) - 1.0);
        cc_ = (4.0 + mu_eff_ / n) / (n + 4.0 + 2.0 * mu_eff_ / n);
        // cc_ = 4.0/(n + 4.0);

        // 确保学习率不会过小
        c1_ = std::max(1e-6, 2.0 / ((n + 1.3) * (n + 1.3) + mu_eff_));
        
        double alpha_mu = 2.0;
        cmu_ = std::min(1.0 - c1_, 
                        alpha_mu * (mu_eff_ - 2.0 + 1.0 / mu_eff_) / ((n + 2.0) * (n + 2.0) + mu_eff_));
        
        // 确保学习率非负
        c1_ = std::max(c1_, 1e-6);
        cmu_ = std::max(cmu_, 1e-6);
    }
    
    Eigen::VectorXd setupWeights(int lambda) 
    {
        int mu = lambda / 2;
        if (mu <= 0) mu = 1;
        
        Eigen::VectorXd weights(mu);
        
        // 标准的对数权重设置
        for (int i = 0; i < mu; ++i) 
        {
            weights(i) = std::log(mu + 0.5) - std::log(i + 1.0);
        }
        
        double sum_weights = weights.sum();
        if (sum_weights > 0) 
        {
            weights /= sum_weights;
        }
        
        return weights;
    }
    
    void updateEigenDecomposition() 
    {
        if (!eigen_decomp_updated_) 
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(covariance_);
            if (es.info() == Eigen::Success) 
            {
                B_ = es.eigenvectors();
                D_ = es.eigenvalues().cwiseSqrt();
                
                // 确保数值稳定性
                double min_eigenvalue = 1e-12;
                double max_eigenvalue = 1e12;
                for (int i = 0; i < D_.size(); ++i) 
                {
                    if (D_(i) < min_eigenvalue) D_(i) = min_eigenvalue;
                    if (D_(i) > max_eigenvalue) D_(i) = max_eigenvalue;
                }
                eigen_decomp_updated_ = true;
            }
        }
    }
    
    void generatePopulation
    (
        std::vector<Eigen::VectorXd>& population,
        std::vector<double>& fitness,
        ThermalPropertyObjective* objective, int n
    ) 
    {
        // 确保特征分解是最新的
        updateEigenDecomposition();
        
        std::normal_distribution<double> normal(0, 1);
        
        for (int i = 0; i < population_size_; ++i) 
        {
            // 生成标准正态随机向量
            Eigen::VectorXd z(n);
            for (int j = 0; j < n; ++j) 
            {
                z(j) = normal(generator_);
            }
            
            // 正确的个体生成公式: x = mean + σ * B * D * z 
            Eigen::VectorXd individual = mean_ + sigma_ * (B_ * D_.asDiagonal() * z);
            
            // 应用边界约束
            std::vector<double> ind_vec(individual.data(), individual.data() + n);
            ind_vec = objective->applyBounds(ind_vec);
            individual = Eigen::Map<Eigen::VectorXd>(ind_vec.data(), n);
            
            population.push_back(individual);
            
            // 计算适应度
            fitness.push_back(objective->evaluate(ind_vec));
        }
    }
    
    void updateDistribution
    (
        const std::vector<Eigen::VectorXd>& population,
        const std::vector<double>& fitness,
        const Eigen::VectorXd& weights
    ) 
    {
        int n = mean_.size();
        int mu = weights.size();
        
        // 对种群按适应度排序
        std::vector<size_t> indices(population_size_);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), 
                 [&](size_t a, size_t b) { return fitness[a] < fitness[b]; });
        
        // 保存旧均值
        Eigen::VectorXd old_mean = mean_;
        
        // 更新均值（加权重组）
        mean_.setZero();
        for (int i = 0; i < mu; ++i) 
        {
            mean_ += weights(i) * population[indices[i]];
        }
        
        // 更新进化路径
        Eigen::VectorXd y_w = (mean_ - old_mean) / sigma_;
        
        // 步长控制路径更新
        Eigen::MatrixXd inv_sqrt_C = B_ * D_.cwiseInverse().asDiagonal() * B_.transpose();
        p_sigma_ = (1.0 - c_sigma_) * p_sigma_ + 
                   std::sqrt(c_sigma_ * (2.0 - c_sigma_) * mu_eff_) * inv_sqrt_C * y_w;
        
        // 判断启发式
        double expectation_norm = std::sqrt(n) * (1.0 - 1.0/(4.0*n) + 1.0/(21.0*n*n));
        double h_sigma_threshold = (1.4 + 2.0/(n + 1.0)) * expectation_norm;
        
        double h_sigma = (p_sigma_.norm() / std::sqrt(1.0 - std::pow(1.0 - c_sigma_, 2.0 * (evaluations_since_restart_ + 1))) 
                         < h_sigma_threshold) ? 1.0 : 0.0;
        
        // 协方差矩阵路径更新
        pc_ = (1.0 - cc_) * pc_ + 
               h_sigma * std::sqrt(cc_ * (2.0 - cc_) * mu_eff_) * y_w;
        
        // 秩1更新
        Eigen::MatrixXd rank1_update = pc_ * pc_.transpose();
        
        // 秩μ更新
        Eigen::MatrixXd rank_mu_update = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < mu; ++i) 
        {
            Eigen::VectorXd y_i = (population[indices[i]] - old_mean) / sigma_;
            rank_mu_update += weights(i) * y_i * y_i.transpose();
        }
        
        // 更新协方差矩阵
        double delta = (1.0 - h_sigma) * cc_ * (2.0 - cc_);
        covariance_ = (1.0 - c1_ - cmu_) * covariance_ + 
                      c1_ * (rank1_update + delta * covariance_) + 
                      cmu_ * rank_mu_update;
        
        // 标记特征分解需要更新
        eigen_decomp_updated_ = false;
        
        // 更新步长
        double norm_p_sigma = p_sigma_.norm();
        double exponent = (c_sigma_ / d_sigma_) * (norm_p_sigma / expectation_norm - 1.0);
        
        // 限制指数范围防止数值问题
        exponent = std::max(-50.0, std::min(5.0, exponent));
        sigma_ *= std::exp(exponent);
        
        // 合理的步长限制
        sigma_ = std::max(1e-12, std::min(1e6, sigma_));
    }
    
    bool checkConvergence(int generation, double current_fitness) 
    {
        // 基础收敛检查：适应度变化很小
        static double last_fitness = 1e100;
        
        if (generation > 200) 
        {
            if (std::abs(last_fitness - current_fitness) < 0.01*config.tolerance) 
            {
                return true;
            }
        }
        last_fitness = current_fitness;
        
        // 检查步长是否过小
        if (sigma_ < 1e-12) 
        {
            return true;
        }
        
        return false;
    }
    
    bool shouldEarlyRestart(int generation, double current_fitness) 
    {
        // 如果长时间没有改善，考虑提前重启
        if (generation > 200 && no_improvement_count_ > 100) 
        {
            return true;
        }
        
        // 如果步长变得非常小但适应度仍然不好，重启
        if (sigma_ < 1e-6 && current_fitness > best_fitness_ * 2) 
        {
            return true;
        }
        
        return false;
    }
    
    int checkGlobalConvergence() 
    {
        // 如果已经达到最佳可能精度
        if (best_fitness_ < config.tolerance) 
        {
            return 1;
        }
        
        // 如果多次重启都没有进步
        if (current_restart_ > 0 && no_improvement_count_ > 0.5*config.max_iterations) 
        {
            return -1;
        }
        
        return 0;
    }
    
public:
    std::string getName() const override 
    { 
        return "协方差矩阵适应进化策略(CMA-ES)"; 
    }

    virtual void outputConfig() const override 
    {
        std::cout << "CMA-ES算法配置：" << std::endl;
        std::cout << "  收敛残差： " << config.tolerance << std::endl;
        std::cout << "  最大重启次数: " << max_restarts_ << std::endl;
        std::cout << "  最大迭代次数: " << config.max_iterations << std::endl;
        std::cout << "  初始种群大小: " << initial_population_size_ << std::endl;
        std::cout << "  初始步长: " << initial_sigma_ << std::endl;
    }

    virtual void setConfig() override 
    {
        std::cout << "输入CMA-ES参数" << std::endl;
        std::cout << "输入最大重启次数: ";
        std::cin >> max_restarts_;
        std::cout << "输入最大迭代次数: ";
        std::cin >> config.max_iterations;
        std::cout << "输入初始种群大小 (0=自动计算): ";
        std::cin >> population_size_;
        initial_population_size_ = population_size_;
        std::cout << "输入初始步长: ";
        std::cin >> sigma_;
        initial_sigma_ = sigma_;
    }
};

// 支持的算法类型
namespace AlgorithmTypes 
{
    const std::string LEVENBERG_MARQUARDT = "LEVENBERG_MARQUARDT";
    const std::string DOGLEG = "DOGLEG";
    const std::string BFGS = "BFGS";
    const std::string LBFGS = "LBFGS";
    const std::string NONLINEAR_CONJUGATE_GRADIENT = "NONLINEAR_CONJUGATE_GRADIENT";
}

// 支持的线性求解器
namespace LinearSolvers 
{
    const std::string DENSE_QR = "DENSE_QR";
    const std::string DENSE_NORMAL_CHOLESKY = "DENSE_NORMAL_CHOLESKY";
    const std::string SPARSE_NORMAL_CHOLESKY = "SPARSE_NORMAL_CHOLESKY";
    const std::string CGNR = "CGNR";
}

// Ceres优化算法
class CeresInversionAlgorithm 
: public InversionAlgorithm 
{
private:
    ceres::Solver::Options solver_options_;
    
public:
    CeresInversionAlgorithm() 
    {
        // 配置求解器选项
        configureSolverOptions();
    }
    
    InversionResult optimize
    (
        ThermalPropertyObjective* objective,
        const std::vector<double>& initial_guess
    ) override 
    {
        InversionResult result;
        config.population_size = initial_guess.size();

        try 
        {
            // 创建Ceres问题
            ceres::Problem problem;
            
            // 使用动态自动微分代价函数
            std::cout << "创建Ceres代价函数..." << std::endl;
            auto* cost_function = new ThermalPropertyCostFunction(objective);
            
            std::vector<std::vector<double>> parameters(2);
            
            // 添加参数块
            size_t k_count = objective->kappa().getValues().size();
            for(size_t i = 0; i < k_count; ++i)
            {
                parameters[0].push_back(initial_guess[i]);
            }
            problem.AddParameterBlock
            (
                parameters[0].data(),
                parameters[0].size()
            );

            // 添加边界约束（归一化热导率）
            for (size_t i = 0; i < k_count; ++i) 
            {
                problem.SetParameterLowerBound(parameters[0].data(), i, 0.0);
                problem.SetParameterUpperBound(parameters[0].data(), i, 1.0);
            }
            
            if(objective->getInversionMode() == 0)
            {
                for(size_t i = 0; i < objective->cp().getValues().size(); ++i)
                {
                    parameters[1].push_back(initial_guess[k_count + i]);
                }

                problem.AddParameterBlock
                (
                    parameters[1].data(),
                    parameters[1].size()
                );

                // 添加边界约束(归一化比热)
                for (size_t i = 0; i < objective->cp().getValues().size(); ++i) 
                {
                    problem.SetParameterLowerBound(parameters[1].data(), i, 0.0);
                    problem.SetParameterUpperBound(parameters[1].data(), i, 1.0);
                }

                // 添加残差块
                std::cout << "添加残差块到Ceres问题..." << std::endl;
                problem.AddResidualBlock
                (
                    cost_function,
                    nullptr,  // 不使用损失函数
                    parameters[0].data(),
                    parameters[1].data()
                );
            }
            else
            {
                // 添加残差块
                std::cout << "添加残差块到Ceres问题..." << std::endl;
                problem.AddResidualBlock
                (
                    cost_function,
                    nullptr,  // 不使用损失函数
                    parameters[0].data()
                );
            }
                    
            // ========== 添加光滑性正则化 ==========
            addSmoothnessRegularization(problem, parameters, objective);
            // =====================================
            
            // 配置线性求解器
            configureLinearSolver(solver_options_);

            // 运行优化
            std::cout<< "运行Ceres求解器..." << std::endl;
            ceres::Solver::Summary summary;
            ceres::Solve(solver_options_, &problem, &summary);
            
            // 处理结果
            result = processSolverResult(summary, parameters, objective);
        } 
        catch (const std::exception& e) 
        {
            result.converged = false;
            result.message = std::string("优化失败: ") + e.what();
        }
        
        return result;
    }
    
    std::string getName() const override 
    {
        return "Ceres优化算法";
    }
    
    virtual void outputConfig() const override
    {
        std::cout << "Ceres算法配置：" << std::endl;
        std::cout << "  最大迭代次数: " << solver_options_.max_num_iterations << std::endl;
        std::cout << "  函数容差: " << solver_options_.function_tolerance << std::endl;
        std::cout << "  梯度容差: " << solver_options_.gradient_tolerance << std::endl;
        std::cout << "  参数容差: " << solver_options_.parameter_tolerance << std::endl;
    }
    
    virtual void setConfig() override
    {
        std::cout << "输入Ceres参数" << std::endl;
        std::cout << "输入最大迭代次数: ";
        std::cin >> solver_options_.max_num_iterations;
        std::cout << "输入函数容差: ";
        std::cin >> solver_options_.function_tolerance;
        std::cout << "输入梯度容差: ";
        std::cin >> solver_options_.gradient_tolerance;
        std::cout << "输入参数容差: ";
        std::cin >> solver_options_.parameter_tolerance;
        std::cout << "选择求解算法：" << std::endl;
        std::cout << "1. LEVENBERG_MARQUARDT\n"
                    "2. DOGLEG\n"
                    "3. BFGS\n"
                    "4. LBFGS\n"
                    "5. NONLINEAR_CONJUGATE_GRADIENT" << std::endl;
        int alg_choice;
        std::cin >> alg_choice;
        if(alg_choice == 1)
        {
            config.algorithm_type = AlgorithmTypes::LEVENBERG_MARQUARDT;
            std::cout<< "选择了LEVENBERG_MARQUARDT算法"<<std::endl;
            std::cout<< "注意：LEVENBERG_MARQUARDT算法适用于小规模问题，参数维度较大时建议选择其他算法以提高效率。" <<std::endl;
            std::cout<< "输入算法参数：" << std::endl;
            std::cout<< "使用单调步长：Y/N " << std::endl;
            char flag;
            std::cin >> flag;
            if(flag == 'Y' || flag == 'y')
            {
                solver_options_.use_nonmonotonic_steps = true;
                std::cout<<"最大非单调步长数：";
                std::cin >> solver_options_.max_consecutive_nonmonotonic_steps;
            }
            std::cout<<"初始置信域半径：";
            std::cin >> solver_options_.initial_trust_region_radius;
            std::cout<<"最大置信域半径：";
            std::cin >> solver_options_.max_trust_region_radius;
        }
    }

private:
    
    // Ceres代价函数封装
    class ThermalPropertyCostFunction 
    : public ceres::CostFunction 
    {
    private:
        ThermalPropertyObjective* objective_;
    public:
        explicit ThermalPropertyCostFunction(ThermalPropertyObjective* objective) 
            : objective_(objective)
        {
            mutable_parameter_block_sizes()->push_back(objective_->kappa().getValues().size());
            if(objective_->getInversionMode() == 0)
            {
                mutable_parameter_block_sizes()->push_back(objective_->extinction().getValues().size());
            }
            
            // 设置残差大小
            set_num_residuals(objective->getResidualCount());
        }
        
        virtual bool Evaluate
        (
            double const* const* parameters,
            double* residuals,
            double** jacobians
        ) const override 
        {
            // 如果需要雅可比矩阵，使用自动微分计算雅可比
            if (jacobians != nullptr) 
            {
                if(jacobians[0] != nullptr)
                {
                    double k_upper = objective_->getThermalConductivityUpperBound();
                    double k_lower = objective_->getThermalConductivityLowerBound();
                    computeJacobianAD(parameters, jacobians[0], k_upper, k_lower, 0);
                }
                
                if(objective_->getInversionMode() == 0 && jacobians[1] != nullptr)
                {
                    double beta_upper = objective_->getExtinctionUpperBound();
                    double beta_lower = objective_->getExtinctionLowerBound();
                    computeJacobianAD(parameters, jacobians[1], beta_upper, beta_lower, 1);
                }
            }

            // 调用目标函数计算残差
            if (residuals != nullptr) 
            {
                objective_->computeResiduals(parameters, residuals);
            }
            
            return true;
        }

    private:
                
        void computeJacobianAD
        (
            const double* const* parameters,
            double* jacobian,
            double param_upper,
            double param_lower, 
            int j
        ) const 
        {
            int num_res = num_residuals();
            int num_param = parameter_block_sizes()[j];
            // 使用Jet类型进行自动微分
            std::vector<ceres::Jet<double, ceres::DYNAMIC>> jet_parameters(num_param);

            // 初始化Jet参数
            for (int i = 0; i < num_param; ++i) 
            {
                jet_parameters[i].a = parameters[j][i];
                jet_parameters[i].v.resize(num_param);
                jet_parameters[i].v.setZero();
                jet_parameters[i].v[i] = 1.0;
            }
            
            // 使用Jet类型计算残差
            std::vector<ceres::Jet<double, ceres::DYNAMIC>> jet_residuals(num_res);
            
            // 初始化Jet残差
            for(int i = 0; i < num_res; ++i)
            {
                jet_residuals[i].a = 0.0;
                jet_residuals[i].v.resize(num_param);
                jet_residuals[i].v.setZero();
            }

            auto ad_solver = createADSolver<ceres::Jet<double, ceres::DYNAMIC>>();

            // 初始化AD求解器
            ad_solver->initialize(num_param);

            // 将归一化参数映射回实际范围
            for (int i = 0; i < num_param; ++i) 
            {
                jet_parameters[i] = jet_parameters[i]*(param_upper - param_lower) + param_lower;
            }

            ad_solver->setupSolverParameters(jet_parameters, j);

            ad_solver->computeResiduals(jet_residuals);

            // 提取雅可比矩阵（行主序存储）
            for (int i = 0; i < num_res; ++i) 
            {
                for (int j = 0; j < num_param; ++j) 
                {
                    jacobian[i * num_param + j] = jet_residuals[i].v[j];
                }
            }
        }

        template<typename T>
        std::unique_ptr<heatTransferSolverAD<T>> 
        createADSolver() const 
        {
            auto& basic_solver = objective_->solver();
            std::string type = basic_solver.typeName();
            std::unique_ptr<heatTransferSolverAD<T>> ad_solver;
            if (type == "singleLayerSolver1D")
            {   
                auto* solverPtr = dynamic_cast<singleLayerSolver1D*>(&basic_solver);
                if (!solverPtr) 
                {
                    throw std::runtime_error("Failed to cast to singleLayerSolver1D");
                }
                ad_solver = std::make_unique<singleLayerSolver1DAD<T>>(*solverPtr);
            }
            else if(type == "multiLayerSolver1D") 
            {
                auto* solverPtr = dynamic_cast<multiLayerSolver1D*>(&basic_solver);
                if (!solverPtr) 
                {
                    throw std::runtime_error("Failed to cast to multiLayerSolver1D");
                }
                ad_solver = std::make_unique<multiLayerSolver1DAD<T>>(*solverPtr);
            } 
            else 
            {
                throw std::runtime_error("不支持的求解器类型");
            }

            return ad_solver;
        }
    };

    // 光滑性正则化代价函数
    class SmoothnessRegularizationCostFunction 
    : public ceres::CostFunction 
    {
    private:
        int num_params_;
        double weight_;
        int order_;  // 1: 一阶差分, 2: 二阶差分
        
    public:
        SmoothnessRegularizationCostFunction(int num_params, double weight = 0.001, int order = 2)
            : num_params_(num_params), weight_(weight), order_(order) 
        {
            // 设置参数块大小
            mutable_parameter_block_sizes()->push_back(num_params);
            
            // 设置残差数量：一阶差分有n-1个残差，二阶差分有n-2个
            if (order_ == 1) 
            {
                set_num_residuals(num_params - 1);
            } 
            else 
            {
                set_num_residuals(num_params - 2);
            }
        }
        
        virtual bool Evaluate
        (
            double const* const* parameters,
            double* residuals,
            double** jacobians
        ) const override 
        {
            const double* params = parameters[0];
            
            if (order_ == 1) 
            {
                // 一阶差分光滑性：最小化相邻参数之差
                for (int i = 0; i < num_params_ - 1; ++i) 
                {
                    residuals[i] = weight_ * (params[i+1] - params[i]);
                }
            } 
            else 
            {
                // 二阶差分光滑性：最小化曲率（更光滑）
                for (int i = 1; i < num_params_ - 1; ++i) 
                {
                    residuals[i-1] = weight_ * (params[i-1] - 2.0 * params[i] + params[i+1]);
                }
            }
            
            // 自动微分计算雅可比矩阵
            if (jacobians != nullptr && jacobians[0] != nullptr) 
            {
                computeSmoothnessJacobian(parameters, jacobians[0]);
            }
            
            return true;
        }
        
    private:
        void computeSmoothnessJacobian
        (
            const double* const* parameters, 
            double* jacobian
        ) const 
        {
            int num_res = num_residuals();
            int num_param = num_params_;
            
            // 使用自动微分计算雅可比
            std::vector<ceres::Jet<double, ceres::DYNAMIC>> jet_params(num_param);
            
            for (int i = 0; i < num_param; ++i) 
            {
                jet_params[i].a = parameters[0][i];
                jet_params[i].v.resize(num_param);
                jet_params[i].v.setZero();
                jet_params[i].v[i] = 1.0;
            }
            
            std::vector<ceres::Jet<double, ceres::DYNAMIC>> jet_residuals(num_res);
            
            if (order_ == 1) 
            {
                for (int i = 0; i < num_res; ++i) 
                {
                    jet_residuals[i] = weight_ * (jet_params[i+1] - jet_params[i]);
                }
            } else 
            {
                for (int i = 0; i < num_res; ++i) 
                {
                    jet_residuals[i] = weight_ * (jet_params[i] - 2.0 * jet_params[i+1] + jet_params[i+2]);
                }
            }
            
            // 提取雅可比矩阵
            for (int i = 0; i < num_res; ++i) 
            {
                for (int j = 0; j < num_param; ++j) 
                {
                    jacobian[i * num_param + j] = jet_residuals[i].v[j];
                }
            }
        }
    };

    void addSmoothnessRegularization
    (
        ceres::Problem& problem,
        std::vector<std::vector<double>>& parameters,
        ThermalPropertyObjective* objective
    ) 
    {
        // 获取配置参数
        double smoothness_weight = objective->getRegularizationParameter();
        bool enable_smoothness = smoothness_weight > 0.0;
        
        if (!enable_smoothness) 
        {
            std::cout << "光滑性正则化已禁用" << std::endl;
            return;
        }
        
        std::cout << "添加光滑性正则化，权重: " << smoothness_weight 
                  << ", 阶数: 2"<< std::endl;
        
        // 热导率参数的光滑性正则化
        int kappa_size = objective->kappa().getValues().size();
        if (kappa_size >= 3) 
        {  
            // 至少需要3个点才能应用二阶差分
            auto* kappa_smooth_cost = 
            new SmoothnessRegularizationCostFunction(kappa_size, smoothness_weight);
            
            problem.AddResidualBlock(kappa_smooth_cost, nullptr, parameters[0].data());

            std::cout << "添加热导率光滑性正则化，参数点数: " << kappa_size << std::endl;
        }
        
        // 衰减参数的光滑性正则化（如果同时反演衰减系数）
        if (objective->getInversionMode() == 0)
        {
            int beta_size = objective->extinction().getValues().size();
            if (beta_size >= 3) 
            {
                auto* beta_smooth_cost = 
                new SmoothnessRegularizationCostFunction(beta_size, smoothness_weight);
                
                problem.AddResidualBlock(beta_smooth_cost, nullptr, parameters[1].data());
                std::cout << "添加比热光滑性正则化，参数点数: " << beta_size << std::endl;
            }
        }
    }

    void configureSolverOptions() 
    {
        ceres::Solver::Options& options = solver_options_;
        
        // 基本配置
        options.max_num_iterations = config.max_iterations;
        options.function_tolerance = 1e-6;
        options.gradient_tolerance = 1e-10;
        options.parameter_tolerance = 1e-6;
        
        // 根据算法类型配置
        if (config.algorithm_type == AlgorithmTypes::LEVENBERG_MARQUARDT) 
        {
            options.minimizer_type = ceres::TRUST_REGION;
            options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
            options.use_nonmonotonic_steps = false;
            options.max_consecutive_nonmonotonic_steps = 5;
            options.initial_trust_region_radius = 1e3;
            options.max_trust_region_radius = 1e16;
        } 
        else if (config.algorithm_type == AlgorithmTypes::DOGLEG) 
        {
            options.minimizer_type = ceres::TRUST_REGION;
            options.trust_region_strategy_type = ceres::DOGLEG;
            options.dogleg_type = ceres::TRADITIONAL_DOGLEG;
        } 
        else if (config.algorithm_type == AlgorithmTypes::BFGS) 
        {
            options.minimizer_type = ceres::LINE_SEARCH;
            options.line_search_direction_type = ceres::BFGS;
            options.line_search_type = ceres::WOLFE;
        } 
        else if (config.algorithm_type == AlgorithmTypes::LBFGS) 
        {
            options.minimizer_type = ceres::LINE_SEARCH;
            options.line_search_direction_type = ceres::LBFGS;
            options.line_search_type = ceres::WOLFE;
        } 
        else if (config.algorithm_type == AlgorithmTypes::NONLINEAR_CONJUGATE_GRADIENT) 
        {
            options.minimizer_type = ceres::LINE_SEARCH;
            options.line_search_direction_type = ceres::NONLINEAR_CONJUGATE_GRADIENT;
            options.line_search_type = ceres::WOLFE;
            options.nonlinear_conjugate_gradient_type = ceres::FLETCHER_REEVES;
        }
    }

    void configureLinearSolver(ceres::Solver::Options& options) 
    {
        if (config.linear_solver == LinearSolvers::DENSE_QR) 
        {
            options.linear_solver_type = ceres::DENSE_QR;
        } 
        else if (config.linear_solver == LinearSolvers::DENSE_NORMAL_CHOLESKY) 
        {
            options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
        } 
        else if (config.linear_solver == LinearSolvers::SPARSE_NORMAL_CHOLESKY) 
        {
            options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        } 
        else if (config.linear_solver == LinearSolvers::CGNR) 
        {
            options.linear_solver_type = ceres::CGNR;
        }
        
        // 根据问题规模智能选择（如果未明确指定）
        if (config.linear_solver.empty()) 
        {
            if (config.population_size > 1000) 
            {
                options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
            } 
            else if (config.population_size > 100) 
            {
                options.linear_solver_type = ceres::DENSE_QR;
            } 
            else 
            {
                options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
            }
        }
    }
    
    InversionResult processSolverResult
    (
        const ceres::Solver::Summary& summary,
        const std::vector<std::vector<double>>& parameters,
        ThermalPropertyObjective* /*objective*/
    ) 
    {
        InversionResult result;
        for(size_t i = 0; i < parameters.size(); ++i)
        {
            for(size_t j = 0; j < parameters[i].size(); ++j)
            {
                result.optimal_parameters.push_back(parameters[i][j]);
            }
        }

        result.optimal_value = summary.final_cost;
        result.iterations = summary.iterations.size();
        result.converged = summary.termination_type == ceres::CONVERGENCE;
        result.message = summary.FullReport();
        
        // 提取成本历史
        for (const auto& iteration : summary.iterations) 
        {
            result.history.push_back(iteration.cost);
        }
        
        return result;
    }
};

// 反演算法工厂
class InversionAlgorithmFactory
{
private:
    std::unordered_map<std::string, std::function<std::unique_ptr<InversionAlgorithm>()>> creators_;
    std::vector<std::string> algorithm_order_; // 保持注册顺序
    
public:
    InversionAlgorithmFactory() 
    {
        // 按热导率反演推荐度排序（索引0为最推荐）
        
        // 索引0: - Ceres优化
        registerAlgorithm
        (
            "Ceres_optimization", 
            []() -> std::unique_ptr<InversionAlgorithm>
            {
                return std::make_unique<CeresInversionAlgorithm>(); 
            }
        );
        
        // 索引1: 最推荐  CMA-ES
        registerAlgorithm
        (
            "cma_es", 
            []() -> std::unique_ptr<InversionAlgorithm> 
            {
                return std::make_unique<CMAES>();
            }
        );
        
        // 索引2: 粒子群优化
        registerAlgorithm
        (
            "particle_swarm", 
            []() -> std::unique_ptr<InversionAlgorithm> 
            {
                return std::make_unique<ParticleSwarmOptimization>();
            }
        );
        
        // 索引3: 遗传算法（不推荐）
        registerAlgorithm
        (
            "genetic_algorithm", 
            []() -> std::unique_ptr<InversionAlgorithm> 
            {
                return std::make_unique<GeneticAlgorithm>();
            }
        );
    }
    
    void registerAlgorithm
    (
        const std::string& name,
        std::function<std::unique_ptr<InversionAlgorithm>()> creator
    ) 
    {
        creators_[name] = creator;
        algorithm_order_.push_back(name); // 记录注册顺序
    }
    
    std::unique_ptr<InversionAlgorithm> createAlgorithm(const std::string& name) 
    {
        auto it = creators_.find(name);
        if (it != creators_.end()) 
        {
            return it->second();
        }
        throw std::invalid_argument("未知算法: " + name);
    }
    
    std::vector<std::string> getAvailableAlgorithms() const 
    {
        return algorithm_order_; // 返回有序列表
    }

    // 获取带描述的算法列表
    std::vector<std::pair<std::string, std::string>> getAlgorithmsWithDescription() const 
    {
        std::vector<std::pair<std::string, std::string>> algo_list;
        
        algo_list.push_back({"Ceres_optimization", "基于Ceres的优化算法 - 适合大规模优化"});
        algo_list.push_back({"cma_es", "CMA-ES - 黑盒优化最强算法"});
        algo_list.push_back({"particle_swarm", "粒子群优化 - 适合连续优化"});
        algo_list.push_back({"genetic_algorithm", "遗传算法 - 适合多峰问题"});

        return algo_list;
    }
};

#endif
