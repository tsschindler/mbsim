#ifndef _WERKZEUGMASCHINE_H
#define _WERKZEUGMASCHINE_H

#include "mbsim/dynamic_system_solver.h"
#include "mbsimFlexibleBody/flexible_body/flexible_body_1s_21_ancf.h"
#include "mbsim/rigid_body.h"
#include <string>

class System : public MBSim::DynamicSystemSolver {
  public:
    System(const std::string &projectName);

  private:
    /** flexible ring */
    MBSimFlexibleBody::FlexibleBody1s21ANCF *rod;

    /** vector of balls */
    std::vector<MBSim::RigidBody*> balls;
};

#endif

