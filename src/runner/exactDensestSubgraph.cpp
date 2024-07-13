#include "../tools/getArgs.hpp"
#include "../densestSubgraph/exactFlowAlgorithm.h"

#include <cassert>
#include <string>
#include <iostream>
#include <ctime>
#include <chrono>
using std::string;

int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);
    
    string filePath;
    if(aC->exist("-f")) filePath = aC->get("-f");
    else {
        std::cout << "-f filePath" << std::endl;
        exit(-1);
    }

    string outFilePath;
    if(aC->exist("-o")) outFilePath = aC->get("-o");

    int p = 4;
    if(aC->exist("-p")) p = atoi(aC->get("-p").c_str());

    int q = 4;
    if(aC->exist("-q")) q = atoi(aC->get("-q").c_str());

    double initialDensity = 0.0;
    if(aC->exist("-i")) initialDensity = std::stof(aC->get("-i"));

    std::cout << filePath << ' ' << outFilePath << std::endl;
    std::cout << "p:" << p << std::endl;
    std::cout << "q:" << q << std::endl;
    std::cout << "initialDensity:" << initialDensity << std::endl;

    exactFlow * runner = new exactFlow(filePath, outFilePath, p, q, initialDensity);
    auto t1 = std::chrono::steady_clock::now();
    
    runner->run();

    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "time:" << duration.count() << "ms" << std::endl;

    delete runner;

    return 0;
}