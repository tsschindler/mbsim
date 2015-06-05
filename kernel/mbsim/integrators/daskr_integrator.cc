/* Copyright (C) 2007  Robert Huber

 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA

 *
 * Contact:
 *   rhuber.berlios.de
 *
 */

#include<config.h>
#include<ctime>
#include<fstream>
#include <mbsim/dynamic_system_solver.h>
#include "fortran/fortran_wrapper.h"
#include "daskr_integrator.h"

#ifndef NO_ISO_14882
using namespace std;
#endif

using namespace fmatvec;
using namespace MBSim;

namespace MBSimIntegrator {

  MBSIM_OBJECTFACTORY_REGISTERXMLNAME(DASKRIntegrator, MBSIMINT % "DASKRIntegrator")

  DASKRIntegrator::DASKRIntegrator() :
      KrylovMethod(0), Ydot(1, INIT, 0), maxOrder(5), info(20, INIT, 0), IntSteps(0), RefusedSteps(0), FuncEvals(0), JacEval(0), RootEvals(0) {
    name = "DASKRIntegrator";
  }

//  void DASKRIntegrator::setDAEIndex(int index) {        // Index 2 or Index 2 with Gear Gupta Leimkuhler Stabilisation (21) or 1
//    DAEIndex = index;                                   // Index 1: DASKR is used as ODE Integrator (la are calculated by mbsim)
//    if (!(DAEIndex == 21)) {
//      assert(DAEIndex<3);
//      assert(DAEIndex>0);
//    }
//  }

  void DASKRIntegrator::setMaxOrder(int order_) {
    maxOrder = order_;
    assert(maxOrder < 6);
    assert(maxOrder > 0);
  }

  // SUBROUTINE RES (T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR)
  void DASKRIntegrator::residuum(double* t, double* Y_, double* Ydot_, double* CJ, double* res_, int* ires, double* rPar, int* iPar) {
    int qSize = iPar[0];
    int uSize = iPar[1];
    //  int xSize = *(iPar+2);
    //  int laSize = *(iPar+3);
    int zSize = iPar[2];
    int YSize = iPar[3];
    int DAEIndex = iPar[5];

    Vec Y(YSize, Y_);
    Vec Ydot(YSize, Ydot_);
    Vec res(YSize, res_);
    cout << "Y =" << Y << endl;
    cout << "system q" << system->getq() << endl;
    cout << "system u" << system->getu() << endl;
    cout << "CJ = " << *CJ << endl;

    system->updateStateDependentVariables(*t);
    system->updateg(*t); // for joint, gap distance need also be updated.
    //checkActive(1);   // todo: add when unilateral constraints are considered.
    system->updategd(*t);
    system->updateT(*t);
    system->updateJacobians(*t); // needed for calculating W.
    system->updateh(*t);
    system->updateM(*t);
    system->updateW(*t);
    //updateV(t);  // only for contact
    //updateG(t);  // not needed.

    system->updatexd(*t);  // needed by res(Ix)
    system->updatewb(*t);  // needed by res(Ila)

    double H = rPar[0];
    double rho = rPar[1];
    double s = nrmInf(system->getM(0)); // todo: d_r , k_r
    double scaling = iPar[4];
    if (!scaling) {
      s = 1;
    }

    cout << "M =" << system->getM(0) << endl;
    cout << "s = " << s << endl;

    // input for residual calculation: Y, Ydot, T, M, h, W, xd(gd), gdd(W, Wb)
    Index Iq(0, qSize - 1);
    Index Iu(qSize, qSize + uSize - 1);
    Index Ix(qSize + uSize, zSize - 1);
    Index Ila(zSize, YSize - 1);
    /* res(Iq) = qdot - T(q)*u
     * res(Iu) = M * udot - h - W * la
     * res(Ix) = xdot - f(x), f(x) = xd;
     * res(Ila) = gdd = W^T * udot + Wdot^T * u  index 1
     */
    res(Iq) = Ydot(Iq) - system->getT() * Y(Iu);
    res(Iu) = system->getM(0) * Ydot(Iu) - system->geth(0) * H * H - system->getW(0) * s * (Y(Ila) + rho * system->getg());  // todo getM(0)?
    res(Ix) = Ydot(Ix) - system->getxd();

    if (DAEIndex == 1){
      res(Ila) = trans(system->getW(0)) * Ydot(Iu) + system->getwb(); // index 1
    }else if(DAEIndex == 2){
      res(Ila) = system->getgd(); // index 2
    }else if (DAEIndex == 3){
      res(Ila) = s * system->getg(); // index 3
    }

//    cout << "t=" << *t << "  " << res << endl;

  }

  void DASKRIntegrator::jac(double* t, double* Y_, double* Ydot_, double* PD, double* CJ, double* rPar, int* iPar) {
//    Vec Y(*(iPar+2), Y_);
//    Mat Jacobian(*(iPar+2), *(iPar+2), PD);
//    system->JacF_DAE(*t, Y, Jacobian, *iPar);
//    Jacobian(0,0,*(iPar+1)-1,*(iPar+1)-1) -= *CJ * Mat(*(iPar+1),*(iPar+1),EYE);
  }

  void DASKRIntegrator::rt(int* neq, double* t, double* Y_, double* Ydot_, int* nrt, double *rval_, double* rPar, int* iPar) {
//    Vec Y(*(iPar+2), Y_);
//    Vec rval(*nrt,rval_);
//    system->root_DAE(Y,rval,*t);
//    //cout << "root evaluation t= "<<*t<<endl;
  }

  void DASKRIntegrator::setFlagErrorTest(int Flag) {
    FlagErrorTest = Flag;
    assert(FlagErrorTest >= 0);           // =0   control errors locally on all variables (diff. and algebraic)
    assert(FlagErrorTest < 2);            // =1   exclude algebraic variables from error test
  }

  void DASKRIntegrator::integratorSetups(DynamicSystemSolver& system_) {
    if (scaling) {
      assert(H > 0);
      tEnd = tEnd / H;
      dtPlot = dtPlot / H;
      dt0 = dt0 / H;
    }else{
      H = 1;
      rho = 0.;
    }

    zSize = system_.getzSize();  // zSize = qSize + uSize[0] + xSize;
    qSize = system_.getqSize();
    uSize = system_.getuSize();
    laBiSize = system_.getlaSize(); // todo: check whether need to different laBiSize and laUniSize

//    int laUniSize= system->getSizeUnilateralConstraints();  // todo: unilateral constarints are considered later.
//    YSize = zSize + laBiSize + laUniSize;
    YSize = zSize + laBiSize;
    // todo: add RT function
//    nrt = system->getsvSize();  // root equations size
    nrt = 0;

    Y.resize(YSize, INIT, 0.0);
    Ydot.resize(YSize, INIT, 0.0);

    t = tStart;
    z >> Y(0, zSize - 1);

    if (z0.size())
      z = z0;
    else
      system_.initz(z);

    // todo:  use static formulation initialize Ydot.
    cout << "initialize Z:" << z << endl;
    cout << "initialize Y:" << Y << endl;
    cout << "initialize Ydot:" << Ydot << endl;

    // aTol and rTol can be chosen both scalars or else both arrays of length NEQ.
    assert(aTol.size() == rTol.size());
    if (aTol.size() == 1)
      iTol = 0; // Scalar
    else {
      iTol = 1; // Vector
      assert(aTol.size() == YSize);
    }

    info(0) = 0;                        // initialisation of integrator
    info(1) = iTol;                     // scalar or vector as tolerances
    info(2) = 0;                        // interval-output mode
    if (dtMax)
      info(6) = 1;              // max stepsize
    if (dt0)
      info(7) = 1;              // initial stepsize
    if (maxOrder != 5)
      info(8) = 1;          // max Order is not default 5, but given by the user. (default for info(8)=0 Order 5)
    info(10) = 1;           // calculate consistent initial values: Given Y_d, calculate Y_a and Y'_d,
    info(13) = 1;                        // stop after initial conditions calculation
    info(15) = FlagErrorTest;            // error test (0: include algebraic variables; 1: exclude algebraic variables)
    info(17) = 2;           // option to get extra printing in initial condition calculation:  0 for NO printing, 1 for minimal printing, or 2 for full printing.

    // direct method
    if (!KrylovMethod) {
      if (useExternalJac)
        info(4) = 1;      // info(4)=1: use external function for jac;
      lrWork = 60 + max(maxOrder + 4, 7) * YSize + 3 * nrt + YSize * YSize + FlagErrorTest * YSize;
      liWork = 40 + YSize + ((info(9) == 1 or info(9) == 3) ? 1 : 0) * YSize + ((info(10) == 1 or info(15) == 1) ? 1 : 0) * YSize;
    }
    // Krylov method
    if (KrylovMethod) {
      assert(LENWP < 0);
      info(11) = 1;      // Krylov method is used to solving the linear system
      int maxl = min(5, YSize);
      lrWork = 60 + (maxOrder + 5) * YSize + 3 * nrt + (maxl + 3) * YSize + (maxl + 3) * maxl + 1 + LENWP + FlagErrorTest * YSize;  // TODO: set LENWP when use krylov.
      liWork = 40 + LENIWP + ((info(9) == 1 or info(9) == 3) ? 1 : 0) * YSize + ((info(10) == 1 or info(15) == 1) ? 1 : 0) * YSize;
    }

    int LID = 40 + ((info(9) == 1 or info(9) == 3) ? 1 : 0) * YSize;

    jroot.resize(nrt, INIT, 0);
    rWork.resize(lrWork, INIT, 0.0);
    iWork.resize(liWork, INIT, 0);

    for (int i = 0; i < zSize; i++)
      iWork(i + LID) = 1;           // differential variable: q, u, x
    for (int i = zSize; i < YSize; i++)
      iWork(i + LID) = -1;           // todo: algebraic variable, is la the only algebraic variable?
    if (info(6))
      rWork(1) = dtMax;                               // max stepsize
    if (info(7))
      rWork(2) = dt0;                                 // initial stepsize
    if (info(8))
      iWork(2) = maxOrder;                            // maximum order
    // not used in the current daskr version
//    iWork(44) = system->getqSize() + system->getxSize();        // Dimension of Index 1 Variables
//    iWork(45) = system->getuSize();                             // Dimension of Index 2 Variables
//    iWork(46) = laBiSize+laUniSize;                             // Dimension of Index 3 Variables

    rPar.resize(2);
    iPar.resize(6, INIT, 0);
    rPar(0) = H;
    rPar(1) = rho;
    iPar(0) = qSize;
    iPar(1) = uSize;
    iPar(2) = zSize;
    iPar(3) = YSize;
    iPar(4) = scaling;
    iPar(5) = DAEIndex;

  }

  void DASKRIntegrator::preIntegrate(DynamicSystemSolver& system_) {
    system = &system_;
    DAEIndex = 1; // default index

    integratorSetups(system_);
    if (FlagCoutInfo) {
      cout << "DAE index " << DAEIndex;
      if (useExternalJac)
        cout << " external jacobian." << endl;
      else
        cout << " internally computed jacobian." << endl;
      if (scaling)
        cout << "Scaling technique will be applied." << endl;
    }

    cout.setf(ios::scientific, ios::floatfield);

    if (FlagPlotIntegrator) {
      integPlot.open((name + ".plt").c_str());
      integPlot << "#1 t [s]:" << endl;
      integPlot << "#2 dt [s]:" << endl;
      integPlot << "#3 order :" << endl;
      integPlot << "#4 idid : " << endl;
      integPlot << "#5 calculation time [s]:" << endl;
    }
  }

  int DASKRIntegrator::computeInitialConditions(MBSim::DynamicSystemSolver& system_, bool FlagPlot, bool FlagCoutInfo_) {
//    info(10)= 1;
//    info(17)= 1;
    // calculate initial conditions
    DDASKR(residuum, &YSize, &t, Y(), Ydot(), &tEnd, info(), rTol(), aTol(), &idid, rWork(), &lrWork, iWork(), &liWork, rPar(), iPar(), jac, 0, rt, &nrt, jroot());
    if (idid == 4) {
      if (FlagCoutInfo_)
        cout << "daskr initial condition calculation was successfull." << endl;
      info(10) = 0;                     // switch off the flag which indicates calculation of initial values
      if (FlagPlot)
        system_.plot(z, t);
    }
    if (idid < 0)
      cout << "DASKR compute Initial Conditions: idid = " << idid << endl;
    return idid;
  }

  void DASKRIntegrator::subIntegrate(MBSim::DynamicSystemSolver& system_, double tEnd) {

    s0 = clock();
    tPlot = t + dtPlot;
    if (tPlot > tEnd)
      tPlot = tEnd;

    while (t < tEnd) {
      // according to the setting (info(3) = 0), the DDASKR may integrate past tOUT and interpolate to obtain the result at tOUT.
      //  and tStop is not defined.
      DDASKR(residuum, &YSize, &t, Y(), Ydot(), &tPlot, info(), rTol(), aTol(), &idid, rWork(), &lrWork, iWork(), &liWork, rPar(), iPar(), jac, 0, rt, &nrt, jroot());

      assert(idid != 1 and idid != 2 and idid != 4);

      // The integration was successfully completed by finding one or more roots of R(T,Y,Y') at T, tPlot may not reached yet.

      //call the code again to continue the integration in the direction of TOUT.
      if (idid == 5) {
        //  todo: maybe need to update the status of the system at this time point.
        // use the old tOut to continue the integration from current t.
      }
      else if (idid == 3) {
        if (scaling) {
          double tReal = t * H;
          Y(qSize, qSize + uSize - 1) = Y(qSize, qSize + uSize - 1) / H;
          system_.plot(z, tReal);

          tPlot = min(tEnd, tPlot + dtPlot);

        }
        else {
          system_.plot(z, t);  // todo: whether plot la.
          // define a new TOUT and continue the integration
          tPlot = min(tEnd, tPlot + dtPlot);
        }
      }
      else if (idid < 0) {
        cout << "DDASKR solver fails at t = " << t << ", the error code is " << idid << endl;
        throw 1;
      }

      if (output)
        cout << "   t = " << t << ",\tdt = " << rWork(6) << ",\torder = " << iWork(7) << "\r" << flush;
      if (FlagPlotIntegrator) {
        integPlot << t << " " << rWork(6) << " " << iWork(7) << " " << idid << " " << time << endl;
      }

    }
    IntSteps += iWork(10);  // todo: to check
    RefusedSteps += iWork(13);
    FuncEvals += iWork(11);
    JacEval += iWork(12);
    RootEvals += iWork(35);

    double s1 = clock();
    time = (s1 - s0) / CLOCKS_PER_SEC;
  }

  void DASKRIntegrator::postIntegrate(MBSim::DynamicSystemSolver& system_) {
    cout.unsetf(ios::scientific);
    if (FlagPlotIntegrator)
      integPlot.close();
    if (FlagPlotIntegrationSum) {
      ofstream integSum((name + ".sum").c_str());
      integSum << "Integration time     : " << time << endl;
      integSum << "Integration steps    : " << IntSteps << endl;
      integSum << "Refused steps        : " << RefusedSteps << endl;
      integSum << "Function calls       : " << FuncEvals << endl;
      integSum << "Jacobian evaluations : " << JacEval << endl;
      integSum << "Root function calls  : " << RootEvals << endl;
      integSum.close();
    }
    if (FlagCoutInfo) {
      if (output)
        cout << endl << endl;
      cout << "Summary Integration with DASKR : " << endl;
      cout << "Integration time     : " << time << endl;
      cout << "Integration steps    : " << IntSteps << endl;
      cout << "Refused steps        : " << RefusedSteps << endl;
      cout << "Function calls       : " << FuncEvals << endl;
      cout << "Jacobian evaluations : " << JacEval << endl;
      cout << "Root function calls  : " << RootEvals << endl;
    }
  }

  void DASKRIntegrator::integrate(MBSim::DynamicSystemSolver& system_) {
    debugInit();
    preIntegrate(system_);
    int errorCodeIniCond = computeInitialConditions(system_, true, FlagCoutInfo);
    if (errorCodeIniCond > 0)
      subIntegrate(system_, tEnd);
    postIntegrate(system_);
  }

}
