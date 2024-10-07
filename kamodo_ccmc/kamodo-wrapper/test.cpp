#include <yaml-cpp/yaml.h>

#include <iostream>

int main() {
    try {
        YAML::Node config = YAML::LoadFile("sample.yaml");

        // The outer element is an array
        for(auto dict : config) {
            // The array element is a map containing the Pos and Rectangle keys:
            auto name = dict["Pos"];
            std::cout << "Name: " << name << '\n';

            auto rect = dict["Rectangle"];

            // loop over the positions Rectangle and print them:
            for(auto pos : rect) {
                std::cout << pos["x"].as<double>() << ",\t"
                          << pos["y"].as<double>() << ",\t"
                          << pos["z"].as<double>() << '\n';
            }
        }

    } catch(const YAML::BadFile& e) {
        std::cerr << e.msg << std::endl;
        return 1;
    } catch(const YAML::ParserException& e) {
        std::cerr << e.msg << std::endl;
        return 1;
    }
}
