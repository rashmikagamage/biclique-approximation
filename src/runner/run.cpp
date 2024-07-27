#include <cassert>
#include <chrono>
#include <ctime>
#include <iostream>
#include <string>

#include "../biClique/BCListPlusPlus.h"
#include "../biClique/BCT.h"
#include "../biClique/BK.h"
#include "../biClique/accuracy.h"
#include "../biClique/bcAndPath.h"
#include "../biClique/colorPath.h"
#include "../biClique/colorPathPequalsQ.h"
#include "../biClique/colorPathSpecificPQ.h"
#include "../biClique/edgePivotSpecificPQ.h"
#include "../biClique/fastEdgePivot.h"
#include "../biClique/pivotAndPath.h"
#include "../biClique/pivotAndPathPequalsQ.h"
#include "../biClique/rawEdgePivot.h"
#include "../biClique/turan.h"
#include "../tools/getArgs.hpp"
using std::string;

int main(int argc, char* argv[]) {
    argsController* aC = new argsController(argc, argv);

    // int deb = 0;
    // if(aC->exist("-deb")) deb = atoi(aC->get("-deb").c_str());

    string filePath;
    if (aC->exist("-f"))
        filePath = aC->get("-f");
    else {
        std::cout << "-f filePath" << std::endl;
        exit(-1);
    }

    string outFilePath;
    if (aC->exist("-o")) outFilePath = aC->get("-o");

    int p = 4;
    if (aC->exist("-p")) p = atoi(aC->get("-p").c_str());

    int q = 4;
    if (aC->exist("-q")) q = atoi(aC->get("-q").c_str());

    std::cout << filePath << ' ' << outFilePath << std::endl;

    // double t1 = clock();

    if (aC->exist("-m")) {
        BK* counter = new BK(filePath, outFilePath);
        auto t1 = std::chrono::steady_clock::now();

        counter->exactCountMaximal(p, q);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;

        delete counter;
    } else if (aC->exist("-pm")) {
        rawEdgePivot* counter = new rawEdgePivot(filePath, outFilePath);

        auto t1 = std::chrono::steady_clock::now();

        if (aC->exist("-v5"))
            counter->exactCountMaximalPivotV2();
        else
            counter->exactCountMaximalPivot();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-apm")) {  // approximate pm
        rawEdgePivot* counter = new rawEdgePivot(filePath, outFilePath);
        double r = 0.1;

        if (aC->exist("-r")) r = atof(aC->get("-r").c_str());

        auto t1 = std::chrono::steady_clock::now();

        counter->approximateCountMaximalPivot(r);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }
    // else if(aC->exist("-fpm")) {//no use delete it
    //     fastEdgePivot * counter = new fastEdgePivot(filePath, outFilePath);
    //     auto t1 = std::chrono::steady_clock::now();

    //     counter->exactCountMaximalPivotFast();

    //     auto t2 = std::chrono::steady_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    //     std::cout << "time:" << duration.count() << "ms" << std::endl;
    //     delete counter;
    // }
    else if (aC->exist("-fpmPQ")) {
        edgePivotSpecificPQ* counter = new edgePivotSpecificPQ(filePath, outFilePath, p, q);
        auto t1 = std::chrono::steady_clock::now();

        if (aC->exist("-bar")) {
            uint32_t bar = atoi(aC->get("-bar").c_str());
            counter->sparseCounting(bar);
        } else
            counter->exactCountMaximalPivotFast();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-cp")) {  // pure color path
        int H = 11;
        if (aC->exist("-H")) H = atoi(aC->get("-H").c_str());
        colorPath* counter = new colorPath(filePath, outFilePath, H, p, q);
        auto t1 = std::chrono::steady_clock::now();
        // counter->approximateCounting();
        // counter->testSubgraphSize();
        uint64_t T = 100000;
        if (aC->exist("-t")) T = atoll(aC->get("-t").c_str());
        // counter->approximateCountingAll(T);

        if (aC->exist("-v5"))
            counter->approximateCountingAllVersion5(T);
        else
            counter->approximateCountingAllVersion2(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-cppq")) {  // pure color path, specific pq
        colorPathSpecificPQ* counter = new colorPathSpecificPQ(filePath, outFilePath, p, q);
        auto t1 = std::chrono::steady_clock::now();
        // counter->approximateCounting();
        // counter->testSubgraphSize();
        uint64_t T = 100000;
        if (aC->exist("-t")) T = atoll(aC->get("-t").c_str());
        // counter->approximateCountingAll(T);

        // if(aC->exist("-v5"))
        //     counter->approximateCountingAllVersion5(T);
        // else
        counter->approximateCountingAllVersion2(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-peqq")) {  // p equals q sampling, p=q<H
        int H = 11;
        if (aC->exist("-H")) H = atoi(aC->get("-H").c_str());
        colorPathPequalsQ* counter = new colorPathPequalsQ(filePath, outFilePath, H, p, q);
        auto t1 = std::chrono::steady_clock::now();
        // counter->approximateCounting();
        // counter->testSubgraphSize();
        uint64_t T = 100000;
        if (aC->exist("-t")) T = atoll(aC->get("-t").c_str());
        // counter->approximateCountingAll(T);

        if (aC->exist("-v5"))
            counter->approximateCountingAllVersion5(T);
        else
            counter->approximateCountingAllVersion2(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-pp")) {  // exact+sampling p<H,q<H
        int H = 11;
        if (aC->exist("-H")) H = atoi(aC->get("-H").c_str());

        pivotAndPath* counter = new pivotAndPath(filePath, outFilePath, H);

        auto t1 = std::chrono::steady_clock::now();
        uint64_t T = 100000;
        if (aC->exist("-t")) T = atoll(aC->get("-t").c_str());

        if (aC->exist("-v5"))
            counter->countingV5(T);
        else
            counter->counting(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-pppeqq")) {  // exact+sampling p=q<H
        int H = 10;
        if (aC->exist("-H")) H = atoi(aC->get("-H").c_str());

        pivotAndPathPequalsQ* counter = new pivotAndPathPequalsQ(filePath, outFilePath, H);

        auto t1 = std::chrono::steady_clock::now();
        uint64_t T = 100000;
        if (aC->exist("-t")) T = atoll(aC->get("-t").c_str());

        if (aC->exist("-v5"))
            counter->countingV5(T);
        else
            counter->counting(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-tu")) {
        turan* counter = new turan(filePath, outFilePath, p, q);
        auto t1 = std::chrono::steady_clock::now();
        uint64_t T = 100000;
        if (aC->exist("-t")) T = atoll(aC->get("-t").c_str());
        counter->sample(T);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-bcpath")) {
        // int H = 11;
        // if(aC->exist("-H")) H = atoi(aC->get("-H").c_str());

        uint64_t T = 100000;
        if (aC->exist("-t")) T = atoll(aC->get("-t").c_str());

        double realV = 1.0;
        if (aC->exist("-v")) realV = atof(aC->get("-v").c_str());

        uint32_t bar = 1000;
        if (aC->exist("-bar")) bar = atoi(aC->get("-bar").c_str());

        bcAndPath* counter = new bcAndPath(filePath, outFilePath, p, q);

        auto t1 = std::chrono::steady_clock::now();

        if (aC->exist("-v5"))
            counter->countingV5(T, realV, bar);
        else
            counter->counting(T, realV, bar);

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    } else if (aC->exist("-my")) {
        double e = 0.01;
        if (aC->exist("e")) e = atof(aC->get("e").c_str());
        if (aC->exist("p")) int p = atoi(aC->get("p").c_str());
        if (aC->exist("q")) int q = atoi(aC->get("q").c_str());
        accuracy* counter = new accuracy(filePath, outFilePath, p, q);

        counter->shadowBuilder1(p, q, e);
    } else if (aC->exist("-zstar")) {
        double e = 0.01;
        if (aC->exist("e")) e = atof(aC->get("e").c_str());
        if (aC->exist("p")) int p = atoi(aC->get("p").c_str());
        if (aC->exist("q")) int q = atoi(aC->get("q").c_str());
        accuracy* counter = new accuracy(filePath, outFilePath, p, q);

        counter->shadowBuilderZStar(p, q, e);
    } else if (aC->exist("-bct")) {
        double e = 0.01;

        BCT* counter = new BCT(filePath, outFilePath, p, q);
        counter->buildTree();

    }

    else {  // BCList++
        BCListPlusPlus* counter = new BCListPlusPlus(filePath, outFilePath, p, q);
        auto t1 = std::chrono::steady_clock::now();

        counter->exactCount();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
        delete counter;
    }

    return 0;
}