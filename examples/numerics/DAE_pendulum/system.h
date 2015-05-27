#ifndef _MASSSPRINGROTATING_H
#define _MASSSPRINGROTATING_H

#include "mbsim/dynamic_system_solver.h"
#include <string>

class System : public MBSim::DynamicSystemSolver {
  public:
    System(const std::string &projectName); 
};

#endif

