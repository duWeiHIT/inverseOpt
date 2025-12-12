#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <stdexcept>

class ConfigReader {
private:
    std::map<std::string, std::string> config_data;
    
public:
    bool loadConfig(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "无法打开配置文件: " << filename << std::endl;
            return false;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            // 跳过空行和注释行
            if (line.empty() || line[0] == '#') continue;
            
            // 查找第一个空格或制表符
            size_t pos = line.find_first_of(" \t");
            if (pos != std::string::npos) {
                std::string key = line.substr(0, pos);
                std::string value = line.substr(pos + 1);
                
                // 移除值前后的空格
                size_t start = value.find_first_not_of(" \t");
                size_t end = value.find_last_not_of(" \t");
                if (start != std::string::npos && end != std::string::npos) {
                    value = value.substr(start, end - start + 1);
                    
                    // 移除行末注释
                    size_t comment_pos = value.find('#');
                    if (comment_pos != std::string::npos) {
                        value = value.substr(0, comment_pos);
                        end = value.find_last_not_of(" \t");
                        if (end != std::string::npos) {
                            value = value.substr(0, end + 1);
                        }
                    }
                }
                
                config_data[key] = value;
            }
        }
        file.close();
        return true;
    }
    
    double getDouble(const std::string& key) const {
        auto it = config_data.find(key);
        if (it != config_data.end()) {
            return std::stod(it->second);
        }
        throw std::runtime_error("配置项未找到: " + key);
    }
    
    int getInt(const std::string& key) const {
        auto it = config_data.find(key);
        if (it != config_data.end()) {
            return std::stoi(it->second);
        }
        throw std::runtime_error("配置项未找到: " + key);
    }
    
    std::string getString(const std::string& key) const {
        auto it = config_data.find(key);
        if (it != config_data.end()) {
            std::string value = it->second;
            // 移除引号
            if (value.front() == '"' && value.back() == '"') {
                value = value.substr(1, value.length() - 2);
            }
            return value;
        }
        throw std::runtime_error("配置项未找到: " + key);
    }
    
    bool getBool(const std::string& key) const {
        auto it = config_data.find(key);
        if (it != config_data.end()) {
            std::string value = it->second;
            return (value == "true" || value == "1" || value == "yes");
        }
        throw std::runtime_error("配置项未找到: " + key);
    }

};

#endif