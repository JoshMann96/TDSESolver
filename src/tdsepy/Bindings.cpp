#include "SimulationManager.h"
#include "KineticOperator.h"
#include "Potentials.h"
#include "Measurers.h"
#include "PhysCon.h"

extern "C" {

    void test(int x){
        std::cout << x << std::endl;
    }

    double test2(int x, double y){
        std::cout << x << y << std::endl;
        return x+y;
    }

    int nottest(){
        std::cout << "nottest ran" << std::endl;
        return 0;
    }


}