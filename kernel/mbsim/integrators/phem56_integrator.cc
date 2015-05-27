/* Copyright (C) 2007 Zhan Wang
 *
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
 *   wangzhanrock@gmail.com
 *
 */

#include <config.h>
#include <mbsim/dynamic_system_solver.h>
#include "fortran/fortran_wrapper.h"
#include "phem56_integrator.h"
#include <time.h>
#include <fstream>

#ifndef NO_ISO_14882
using namespace std;
#endif

using namespace fmatvec;
using namespace MBSim;
using namespace MBXMLUtils;
using namespace xercesc;

namespace MBSimIntegrator {

  MBSIM_OBJECTFACTORY_REGISTERXMLNAME(PHEM56Integrator, MBSIMINT % "PHEM56Integrator")

  PHEM56Integrator::PHEM56Integrator() :
      mode(0) {
    name = "PHEM56Integrator";
  }

  double PHEM56Integrator::tPlot = 0;
  double PHEM56Integrator::dtOut = 0;
  Vec PHEM56Integrator::zInp;
  ofstream PHEM56Integrator::integPlot;
  double PHEM56Integrator::s0;
  double PHEM56Integrator::time = 0;
  bool PHEM56Integrator::output_;
  std::vector<fmatvec::Vec2> PHEM56Integrator::la_plot;
  std::ofstream PHEM56Integrator::la0Sum;
  std::ofstream PHEM56Integrator::la1Sum;
  Vec PHEM56Integrator::rPar;
  Vector<Var, int> PHEM56Integrator::iPar;

  /**
   * NBLK: NUMBER OF BLOCKS OF AM
   * NMRC (SIZE OF A BLOCK OF AM)
   * NPGP (0 IF GP AS THE SAME PATTERN AS PREVIOUS CALL)
   * NPFL (0 IF FL AS THE SAME PATTERN AS PREVIOUS CALL)
   */

//  void PHEM56Integrator::fprob(int* IFCN, int* qSize, int* uSize, int* xSize, int* laSize, int* LDG, int* LDF, int* LDA, int* NBLK, int* NMRC, int* NPGP, int* NPFL, int* INDGR, int* INDGC, int* INDFLR, int* INDFLC, double* t, double* q_, double* u_, double* x_, double* la_, double* g, double* GQ, double* h, double* GQQ, double* w_hat, double* FL, double* qd, double* xd, double* M) {
//
//    int m = system->getW(0).cols();
//    int n = system->getW(0).rows();
//
//    int ifcn = (int)(*IFCN);
//    Mat W_mbsim_neg;
////    q = Vec(*qSize, q_);
////    u = Vec(*uSize, u_);
////    x = Vec(*xSize, x_);
////    la = Vec(*laSize, la_);
//    system->setq(Vec(*qSize, q_));
//    system->setu(Vec(*uSize, u_));
//    system->setx(Vec(*xSize, x_));
////    system->setLa(Vec(*laSize, la_));
//
//
//    cout << "-----------case " << *IFCN << "------------------" << endl;
//    cout << "time =" << *t << endl;
//    cout << "q" << Vec(*qSize, q_) << endl;
//    cout << "u" << Vec(*uSize, u_) << endl;
//
//    if ((ifcn == 1) || (ifcn >= 7)) // compute Mass AM
//    {
//      system->updateM(*t);
//      system->getM(0).copy(M);
//    }
//    if ((ifcn == 1) || (ifcn == 5) || (ifcn == 7) || (ifcn == 8)) // compute F
//    {
//      // in order to update h, the state dependent variables, gap distance and Jacobian matrices have to be updated first
//      system->updateStateDependentVariables(*t);
//      system->updateg(*t); // h vector may contains the spring forces which evaluation needs the updated gap distance.
//      system->updateJacobians(*t);
//      system->updateh(*t);
//      system->geth().copy(h);
//
//      cout << "h_mbsim =" << system->geth() << endl;
//      cout << " h =  " << endl;
//      for (int i = 1; i < *uSize; i++)
//        cout << h[i] << endl;
//    }
//    if (ifcn == 4) // compute G (constraints)
//    {
//      system->updateg(*t);
//      Vec g_neg = -system->getg();
//      g_neg.copy(g); // todo???
////      system->getg().copy(g); // todo???
//    }
//
//    if ((ifcn == 1) || (ifcn == 6) || (ifcn >= 10))  // compute GQ ( Jacobian of the constraints)
//    {
//      system->updateJacobians(*t);
//      system->updateW(*t);
//      W_mbsim_neg = -system->getW(0);
//      W_mbsim_neg.T().copy(GQ);
//
//    }
//
//    if ((ifcn == 5) || (ifcn == 7))  // compute GPP ( Hessian of the constraints)
//    {
//      //RuntimeException::selfThrow("Hem5OSI::fprob(), G_qq is not available");
//      *IFCN = 0;
//    }
//
//    if ((ifcn == 3) || (ifcn == 6) || (ifcn >= 10))  // compute GT (partial time derivative of the constraints)
//    {
//
//    }
//
//    if (ifcn == 0) // compute UDOT
//    {
//      system->updatexd(*t);
//      system->getxd().copy(xd);  // copy the ele memeory of MBSim::Vec::xd to the double* xd which is used by fortran subroutine.
//                                 // the memory of the array where double* xd points to locates in WK(IXD).
//    }
//
//    if ((ifcn == 1) || (ifcn == 2) || (ifcn == 10))  // compute QDOT
//    {
//      system->updateT(*t);
//      system->updateqd(*t);
//      system->getqd().copy(qd);
//    }
//  }
  void PHEM56Integrator::fprob(int* IFCN, int* qSize, int* uSize, int* xSize, int* laSize, int* LDG, int* LDF, int* LDA, int* NBLK, int* NMRC, int* NPGP, int* NPFL, int* INDGR, int* INDGC, int* INDFLR, int* INDFLC, double* t, double* q_, double* u_, double* x_, double* la_, double* g, double* GQ, double* h, double* GQQ, double* GT, double* FL, double* qd, double* xd, double* M) {

    int ifcn = (int) (*IFCN);
    Mat W_mbsim_neg;
    //    q = Vec(*qSize, q_);
    //    u = Vec(*uSize, u_);
    //    x = Vec(*xSize, x_);
    //    la = Vec(*laSize, la_);
    system->setq(Vec(*qSize, q_));
    system->setu(Vec(*uSize, u_));
    system->setx(Vec(*xSize, x_));
    //    system->setLa(Vec(*laSize, la_));

    system->updateStateDependentVariables(*t);
    system->updateg(*t); // h vector may contains the spring forces which evaluation needs the updated gap distance.

//    system->updategd(*t); // todo check!

    system->updateT(*t);
    system->updateJacobians(*t);
    system->updateh(*t);
    system->updateM(*t);
    system->updateW(*t);

    system->updatexd(*t);
    system->updateqd(*t);

    int m = system->getW(0).rows();
    int n = system->getW(0).cols();

    cout << "-----------case " << *IFCN << "------------------" << endl;
    cout << "time =" << *t << endl;
    cout << "q" << Vec(*qSize, q_) << endl;
    cout << "u" << Vec(*uSize, u_) << endl;

    if ((ifcn == 1) || (ifcn >= 7)) // compute Mass AM
    {
//        system->updateM(*t);
      system->getM(0).copy(M);
      cout << "M mbsim =" << system->getM(0) << endl;
      cout << "M phem56 = " << endl;
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          cout << M[j * m + i] << " ";
        }
        cout << endl;
      }
    }
    if ((ifcn == 1) || (ifcn == 5) || (ifcn == 7) || (ifcn == 8)) // compute F
    {
      // in order to update h, the state dependent variables, gap distance and Jacobian matrices have to be updated first
//        system->updateStateDependentVariables(*t);
//        system->updateg(*t); // h vector may contains the spring forces which evaluation needs the updated gap distance.
//        system->updateJacobians(*t);
//        system->updateh(*t);
      system->geth().copy(h);

      cout << "h_mbsim =" << system->geth() << endl;
      cout << " h phem56 =  " << endl;
      for (int i = 0; i < *uSize; i++)
        cout << h[i] << endl;
    }
    if (ifcn == 4) // compute G (constraints)
    {
//        system->updateg(*t);
      Vec g_neg = -system->getg();
      g_neg.copy(g); // todo???
//      system->getg().copy(g); // todo???
    }

    if ((ifcn == 1) || (ifcn == 6) || (ifcn >= 10))  // compute GQ ( Jacobian of the constraints)
    {
//        system->updateJacobians(*t);
//        system->updateW(*t);
      W_mbsim_neg = -system->getW(0);
      W_mbsim_neg.T().copy(GQ);

      cout << "W mbsim =" << system->getW(0) << endl;
      cout << "W_mbsim_neg =" << W_mbsim_neg << endl;
      cout << "GQ phem56 = " << endl;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
          cout << GQ[j * n + i] << " ";
        }
        cout << endl;
      }

    }

    if ((ifcn == 5) || (ifcn == 7))  // compute GPP ( Hessian of the constraints)
    {
      //RuntimeException::selfThrow("Hem5OSI::fprob(), G_qq is not available");
      *IFCN = 0;
    }

    if ((ifcn == 3) || (ifcn == 6) || (ifcn >= 10))  // compute GT (partial time derivative of the constraints)
    {
      for (int i = 0; i < *laSize; i++)
        GT[0] = 0.0;
    }

    if (ifcn == 0) // compute UDOT
    {
//        system->updatexd(*t);
      system->getxd().copy(xd);  // copy the ele memeory of MBSim::Vec::xd to the double* xd which is used by fortran subroutine.
                                 // the memory of the array where double* xd points to locates in WK(IXD).
      cout << "xd_mbsim =" << system->getxd() << endl;
      cout << " xd phem56 =  " << endl;
      for (int i = 0; i < *xSize; i++)
        cout << xd[i] << endl;
    }

    if ((ifcn == 1) || (ifcn == 2) || (ifcn == 10))  // compute QDOT
    {
//        system->updateT(*t);
//        system->updateqd(*t);
      system->getqd().copy(qd);
      cout << "qd =" << system->getqd() << endl;
      cout << " qd phem56 =  " << endl;
      for (int i = 0; i < *uSize; i++)
        cout << qd[i] << endl;
    }
  }

  /* fprob1 is uesd to test the directly implemented simple pendulum */
  void PHEM56Integrator::fprob1(int* IFCN, int* qSize, int* uSize, int* xSize, int* laSize, int* LDG, int* LDF, int* LDA, int* NBLK, int* NMRC, int* NPGP, int* NPFL, int* INDGR, int* INDGC, int* INDFLR, int* INDFLC, double* t, double* q_, double* u_, double* x_, double* la_, double* g, double* GQ, double* h, double* GQQ, double* w_hat, double* FL, double* qd, double* xd, double* M) {

    int m = 2;
    int n = 1;

    Vec q = Vec(*qSize, q_);
    Vec u = Vec(*uSize, u_);
    Vec x = Vec(*xSize, x_);
//    Vec la = Vec(*laSize, la_);

//    Vec qd_m;

//    Vec h_m(2, INIT);
//    Mat GQ_m(1, 2);

    cout << "-----------case " << *IFCN << "------------------" << endl;
    cout << "q" << q << endl;
    cout << "u" << u << endl;

    int ifcn = (int) (*IFCN);

    if ((ifcn == 1) || (ifcn >= 7)) // compute Mass AM
    {
      SymMat M_m(2, INIT);
      M_m(0, 0) = 1.0;
      M_m(1, 1) = 1.0;
      M_m.copy(M);

      cout << "M mbsim =" << M_m << endl;
      cout << "M phem56 = " << endl;
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
          cout << M[j * m + i] << " ";
        }
        cout << endl;
      }

    }
    if ((ifcn == 1) || (ifcn == 5) || (ifcn == 7) || (ifcn == 8)) // compute F
    {
      Vec h_m(2, INIT);
      h_m(1) = -9.81;
      h_m.copy(h);

      cout << "h_mbsim =" << system->geth() << endl;
      cout << " h phem56 =  " << endl;
      for (int i = 0; i < *uSize; i++)
        cout << h[i] << endl;
    }
    if (ifcn == 4) // compute G (constraints)
    {
      g[0] = q(0) * q(0) + q(1) * q(1) - 0.01;
      cout << "g = " << g[0];
    }

    if ((ifcn == 1) || (ifcn == 6) || (ifcn >= 10))  // compute GQ ( Jacobian of the constraints)
    {
      Mat GQ_m(1, 2);
      GQ_m = -2. * q.T(); // W_trans = 2*q.T();
      GQ_m.copy(GQ);

      cout << "W mbsim =" << 2. * q << endl;
      cout << "W_mbsim_neg =" << -2. * q << endl;
      cout << "GQ phem56 = " << endl;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
          cout << GQ[j * n + i] << " ";
        }
        cout << endl;
      }

    }

    if ((ifcn == 5) || (ifcn == 7))  // compute GPP ( Hessian of the constraints)
    {
      *IFCN = 0;
    }

    if ((ifcn == 3) || (ifcn == 6) || (ifcn >= 10))  // compute GT (partial time derivative of the constraints)
    {

    }

    if (ifcn == 0) // compute UDOT
    {

    }

    if ((ifcn == 1) || (ifcn == 2) || (ifcn == 10))  // compute QDOT
    {
      Vec qd_m;
      qd_m = u;
      qd_m.copy(qd);

      cout << "qd =" << qd_m << endl;
      cout << " qd phem56 =  " << endl;
      for (int i = 0; i < *uSize; i++)
        cout << qd[i] << endl;
    }
  }

  void PHEM56Integrator::plot(int* NSTEP, int* qSize, int* uSize, int* xSize, int* laSize, int* LRDO, double* q, double* u, double* x, double* ud, double* la, double* DOWK) {

    double tReal, dt;
    double H = rPar(0);

    double t = DOWK[0] + DOWK[1];  // tOld=DOWK[0]; dt=DOWK[1] ;
    dt = DOWK[1];
    while (t >= tPlot) {
      double zSize = *qSize + *uSize + *xSize;
      int FIRST = 1;
      for (int i = 1; i <= zSize; i++)
        zInp(i - 1) = POL4(&i, &FIRST, qSize, uSize, xSize, LRDO, &tPlot, DOWK); // POL4(I,FIRST,NQ,NV,NU,LRDO,X,DOWK)

      tReal = tPlot * H;
      zInp(*qSize, *qSize + *uSize - 1) = zInp(*qSize, *qSize + *uSize - 1) / H;
      dt = dt * H;

      system->plot(zInp, tReal);
      if (output_)
        cout << "   t = " << tReal << ",\tdt = " << dt << "\r" << flush;

//      la_plot.push_back(Vec(*nla, la));
      la0Sum << t << " " << la[0] << endl;
      la1Sum << t << " " << la[1] << endl;

      double s1 = clock();
      time += (s1 - s0) / CLOCKS_PER_SEC;
      s0 = s1;

      integPlot << tReal << " " << dt << " " << time << endl;
      tPlot += dtOut;
    }
  }

  void PHEM56Integrator::integratorSetups(DynamicSystemSolver& system_) {
    if (scaling) {
      assert(H > 0);
      tEnd = tEnd / H;
      dtPlot = dtPlot / H;
      dt0 = dt0 / H;
    }
    else {
      H = 1;
      rho = 0.;
    }

    qSize = system_.getqSize();
    uSize = system_.getuSize();
    xSize = system_.getxSize();
    laBiSize = system_.getlaSize();
    zSize = qSize + uSize + xSize;
    YSize = zSize + laBiSize;

    Index Iq(0, qSize - 1);
    Index Iu(qSize, qSize + uSize - 1);
    Index Ix(qSize + uSize, zSize - 1);
    Index Iz(0, zSize - 1);
    Index Ila(zSize, YSize - 1);
    ud.resize(uSize, INIT, 0.0);
    Y.resize(YSize, INIT, 0.0);
    q >> Y(Iq);
    u >> Y(Iu);
    x >> Y(Ix);
    z >> Y(Iz);
    la >> Y(Ila);

    if (z0.size())
      z = z0;
    else
      system_.initz(z);

    //  // TODO: for the hard implemented simple pendulum example, if it is used, in the subIntegrate function, fprob1 should be used instead of fprob.
//        zSize = 4;
//        qSize = 2;
//        uSize = 2;
//        xSize = 0;
//        laBiSize = 1;
//        YSize = zSize + laBiSize;
//        Index Iq(0, qSize - 1);
//        Index Iu(qSize, qSize + uSize - 1);
//        Index Ix(qSize + uSize, zSize - 1);
//        Index Iz(0, zSize - 1);
//        Index Ila(zSize, YSize - 1);
//        Y.resize(YSize, INIT, 0.0);
//        ud.resize(uSize, INIT, 0.0);
//        q >> Y(Iq);
//        u >> Y(Iu);
//        x >> Y(Ix);
//        z >> Y(Iz);
//        la >> Y(Ila);
//        double L = 0.2;
//        double phi = -M_PI / 6;
//        Y(0) = 0.5 * L * cos(phi);
//        Y(1) = 0.5 * L * sin(phi);
    //      Y(4) = 9.81;
    //      aTol.resize(4);
    //      rTol.resize(4);
    //      aTol(0) = 1e-4;
    //      aTol(1) = 1e-4;
    //      aTol(2) = 1e-5;
    //      aTol(3) = 1e-5;
    //      rTol(0) = 1e-3;
    //      rTol(1) = 1e-3;
    //      rTol(2) = 1e-4;
    //      rTol(3) = 1e-4;

    t = tStart;
    assert(aTol.size() == rTol.size());
    if (aTol.size() == 1) {
      iTol = 0; // Scalar
    }
    else {
      iTol = 1; // Vector
      assert(aTol.size() == zSize);
    }

    iout = 2; // dense output is required

    rPar.resize(2);
    iPar.resize(1);
    rPar(0) = H;
    rPar(1) = rho;
    iPar(0) = scaling;

    int LL, NZA;

    // todo: parameters needed fur sparse mode
    int LDG, LDF, NMRC, NBLK, IXS, IS;
    LDG = 0;
    LDF = 0;
    NMRC = uSize;
    NBLK = 1;
    IXS = 0;
    NZA = 0;
    IS = 0;

    if (mode <= 3) {
      LL = 8 * laBiSize * uSize + 4 * (laBiSize + uSize) * (laBiSize + uSize);  // full matrix mode;
      LDG = laBiSize;
      LDF = laBiSize;
      NZA = LDG + max(LDG, LDF) + NMRC * NMRC * NBLK;
    }
    else if (mode == 4 || mode == 5) {
      LL = 6 * LDG + 2 * LDF + 4 * NMRC * NMRC * NBLK; //matrix matrix mode
      NZA = LDG + max(LDG, LDF) + NMRC * NMRC * NBLK;
    }

    // todo: check the size;

    lrWork = 13 + 16 * qSize + 15 * uSize + 15 * xSize + 4 * laBiSize + 2 * NZA + 2 * IXS + LL;
    liWork = 95 + 2 * (uSize + laBiSize) + 2 * IS + 12 * LDG + 4 * LDF + 4 * NZA;
    //  int lWork = 1200;
    //  int liWork = 1200;

    rWork.resize(lrWork, INIT);
    iWork.resize(liWork, INIT);

    if (dtMax > 0)
      rWork(4) = dtMax;

    iWork(10) = maxSteps;   //Maximum Step Numbers
    iWork(11) = 1; // initial projection //
    iWork(12) = 1; //4; // a projection is done every NBS accepted steps.
    iWork(13) = mode;
    iWork(14) = 1; // at the final stage of every time step, projection onto the index 1 constraints is performed.

  }

  void PHEM56Integrator::preIntegrate(DynamicSystemSolver& system_) {
    system = &system_;
    DAEIndex = 2; // default index

    integratorSetups(system_);
    if (FlagCoutInfo) {
      cout << "DAE index " << DAEIndex;
      if (useExternalJac)
        cout << " external jacobian." << endl;
      else
        cout << " internally computed jacobian." << endl;
      if (scaling)
        cout << "Scaling technique will be applied. Scaling factor is: H = " << H << "." << endl;
    }

    cout.setf(ios::scientific, ios::floatfield);

    if (FlagPlotIntegrator) {
      integPlot.open((name + ".plt").c_str());
      integPlot << "#1 t [s]:" << endl;
      integPlot << "#2 dt [s]:" << endl;
      integPlot << "#3 idid : " << endl;
      integPlot << "#4 calculation time [s]:" << endl;
    }
  }

  void PHEM56Integrator::subIntegrate(MBSim::DynamicSystemSolver& system_, double tEnd) {
    tPlot = t + dtPlot;
    dtOut = dtPlot;
    system_.plot(z, t);

    zInp.resize(zSize);
    output_ = output;

    s0 = clock();

    la0Sum.open("la0.out");
    la1Sum.open("la1.out");

    // the memmory of q, u, x, ud,, la are created in this c++ driver function. Other internal array and matrix (eg. G, M, W_trans and so on) are created internally by the FORTRAN subroutine.
    if (1)
      PHEM56(&qSize, &uSize, &xSize, &laBiSize, fprob, &t, q(), u(), x(), ud(), la(), &tEnd, &dt0, rTol(), aTol(), &iTol, plot, &iout, rWork(), &lrWork, iWork(), &liWork, &idid);
    else
      PHEM56(&qSize, &uSize, &xSize, &laBiSize, fprob1, &t, q(), u(), x(), ud(), la(), &tEnd, &dt0, rTol(), aTol(), &iTol, plot, &iout, rWork(), &lrWork, iWork(), &liWork, &idid);
  }

  void PHEM56Integrator::postIntegrate(MBSim::DynamicSystemSolver& system_) {
    cout.unsetf(ios::scientific);
    if (FlagPlotIntegrator)
      integPlot.close();
    if (FlagPlotIntegrationSum) {
      ofstream integSum((name + ".sum").c_str());
      integSum << "Integration time: " << time << endl;
//          integSum << "Integration steps: " << integrationSteps << endl;
      integSum.close();
    }

    la0Sum.close();
    la1Sum.close();

    cout << endl;

  }

}
