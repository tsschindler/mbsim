/* Copyright (C) 2004-2006  Martin Förg, Robert Huber

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
 *   mfoerg@users.berlios.de
 *
 */

#include<config.h>
#include<ctime>
#include<fstream>
#include <mbsim/dynamic_system_solver.h>
#include "fortran/fortran_wrapper.h"
#include "radau5dae_integrator.h"

#ifndef NO_ISO_14882
using namespace std;
#endif

using namespace MBSim;
using namespace fmatvec;

namespace MBSimIntegrator {

  MBSIM_OBJECTFACTORY_REGISTERXMLNAME(RADAU5DAEIntegrator, MBSIMINT % "RADAU5DAEIntegrator")

  RADAU5DAEIntegrator::RADAU5DAEIntegrator() :
      updateJacobian(0) {
    name = "RADAU5DAEIntegrator";
    DAEIndex = 3; // default index when scaling technique is applied.
  }

  double RADAU5DAEIntegrator::tPlot = 0;
  double RADAU5DAEIntegrator::dtOut = 0;
  Vec RADAU5DAEIntegrator::YInp;
  std::ofstream RADAU5DAEIntegrator::integPlot;
  double RADAU5DAEIntegrator::s0;
  double RADAU5DAEIntegrator::time = 0;
  bool RADAU5DAEIntegrator::output_;

  void RADAU5DAEIntegrator::setDAEIndex(int index) {    // Index 2 or 3 or Index 2 with Gear Gupta Leimkuhler Stabilisation (21)
    DAEIndex = index;
    assert(DAEIndex < 4);
    assert(DAEIndex > 0);
  }

//  void RADAU5DAEIntegrator::fdae(int* YSize, double* t, double* Y_, double* F_, double* rpar, int* ipar) {
//    Vec Y(*YSize, Y_);
//    Vec F(*YSize, F_);
////    system->F_DAE(Y, F, *t, *ipar);     //ipar[0]: DAEIndex
//    int DAEIndex = ipar[0];
//    int zSize = ipar[1];
////    int laBilateralSize = ipar[3];
//
////    int IndexlaUnilaterlal;
//    // originally F(0, zSize -1) is updated via updating zd.
////    if (qd() != F()) {
////      zdParent >> F(0, zSize - 1);
////      updatezdRef (zdParent);
////    }
////    if (q() != YParent()) {
////      updatezRef(YParent(0, zSize - 1));
////    }
//
//    system->updateStateDependentVariables(*t);
//    system->updateg(*t);
//    //checkActive(1);   // todo: add when unilateral constraints are considered.
//    system->updategd(*t);
//    system->updateT(*t);  // updateg, updategd, updateT are from Dynamicsystem.
//    system->updateJacobians(*t); // needed for calculating W.
//    system->updateh(*t);
//    system->updateM(*t);
//    system->updateW(*t);
//    system->facLLM(); // from MBSim::Group
//
//    // todo: setVauledForce law?
////    if (linkSetValuedBilateral.size() || linkSetValuedUnilateral.size()) {
////      updateW(t);
////      updateGb(t);
////      updater(t);
////    }
//    // todo: replace  F(Iq), F(Iu) updatezd by direct formulation?
//    // the M^-1 problem?
//    system->updatezd(*t);  // todo: call the group function, but get function call the dynamic system.
//
//    int qSize = system->getqSize();
//    int uSize = system->getuSize();
//    Index Iq(0, qSize - 1);
//    Index Iu(qSize, qSize + uSize-1);
//    Index Ix(qSize + uSize, zSize-1);
//    Index Ila(zSize, *YSize - 1);
//
//    F(Iq) = system->getqd(); // todo: call the dynamic_system qd;
//    F(Iu) = system->getud();
//    F(Ix) = system->getxd();
//    if (DAEIndex == 2) {
//      F(Ila) = system->getgd();
//    }
//    if (DAEIndex == 3)
//      F(Ila) = system->getg(); //TODO Unilateral
//
//    cout << "t ="<<*t<< endl <<"  F= "<<F<<endl;
//    cout << "Y =" << Y << endl;
//    cout << "M =" << system->getM(0)<<endl;;
//    cout << "h =" << system->geth(0)<<endl;;
//    cout << "W = "<< system->getW(0)<<endl;;
//
////    if (DAEIndex == 21) { // Gear Gupta Leimkuhler Formulation (GGL)        // TODO Unilateral
////      F(zSize, zSize + laBilateralSize - 1) = g(gIndBilateral);
////      F(zSize + laBilateralSize, zSize + laBilateralSize + laBilateralSize - 1) = gd(laIndBilateral);
////      qd += W(Index(0, uSize - 1), laIndBilateral) * YParent(zSize + laBilateralSize, zSize + laBilateralSize + laBilateralSize - 1);
////    }
//  }

  void RADAU5DAEIntegrator::fdae(int* YSize, double* t, double* Y_, double* F_, double* rPar, int* iPar) {
    Vec Y(*YSize, Y_);
    Vec F(*YSize, F_);
    //    system->F_DAE(Y, F, *t, *ipar);     //ipar[0]: DAEIndex
    int DAEIndex = iPar[0];
    int zSize = iPar[1];

    system->updateStateDependentVariables(*t);
    system->updateg(*t);
    system->updategd(*t);
    system->updateT(*t);  // updateg, updategd, updateT are from Dynamicsystem.
    system->updateJacobians(*t); // needed for calculating W.
    system->updateh(*t);
    system->updateM(*t);
    system->facLLM(); // from MBSim::Group

    system->updateW(*t);

    system->updatexd(*t);  // needed by F(Ix)
    system->updatewb(*t);  // needed by F(Ila)

    int qSize = system->getqSize();
    int uSize = system->getuSize();

    double H = rPar[0];
    double rho = rPar[1];
    double s = nrmInf(system->getM(0)); // todo: d_r , k_r
    double scaling = iPar[4];
    if (!scaling) {
      s = 1;
    }

    cout << "M =" << system->getM(0) << endl;
    cout << "s = " << s << endl;

    Index Iq(0, qSize - 1);
    Index Iu(qSize, qSize + uSize - 1);
    Index Ix(qSize + uSize, zSize - 1);
    Index Ila(zSize, *YSize - 1);

    F(Iq) = system->getT() * Y(Iu);
    F(Iu) = slvLLFac(system->getLLM(0), system->geth(0) * H * H + system->getW(0) * s * (Y(Ila) + rho * system->getg()));
    F(Ix) = system->getxd();

//    system->zdot(Y, F, *t);

    if (DAEIndex == 2) {
//      F(Ila) = system->getgd();
      F(Ila) = trans(system->getW(0)) * Y(Iu);
    }
    if (DAEIndex == 3)
      F(Ila) = s * system->getg();
    if (DAEIndex == 1)
      F(Ila) = trans(system->getW(0)) * F(Iu) + system->getwb();

//    cout << "t =" << *t << endl << "  F= " << F << endl;
//    cout << "Y =" << Y << endl;
//    cout << "M =" << system->getM(0) << endl;
//    cout << "h =" << system->geth(0) << endl;
//    cout << "W = " << system->getW(0) << endl;

    //    if (DAEIndex == 21) { // Gear Gupta Leimkuhler Formulation (GGL)        // TODO Unilateral
    //      F(zSize, zSize + laBilateralSize - 1) = g(gIndBilateral);
    //      F(zSize + laBilateralSize, zSize + laBilateralSize + laBilateralSize - 1) = gd(laIndBilateral);
    //      qd += W(Index(0, uSize - 1), laIndBilateral) * YParent(zSize + laBilateralSize, zSize + laBilateralSize + laBilateralSize - 1);
    //    }
  }

  void RADAU5DAEIntegrator::mdae(int* YSize, double* MasMat_, int* lMas, double* rPar, int* iPar) {
    Mat MasMat(*lMas, *YSize, MasMat_);
    int zSize = *(iPar + 1);
    for (int i = 0; i < zSize; i++)
      MasMat(0, i) = 1.0;
    for (int i = zSize; i < *YSize; i++)
      MasMat(0, i) = 0.0;
  }

  void RADAU5DAEIntegrator::jac(int* YSize, double* t, double* Y_, double* jac_, int* ljac_, double* rPar, int* iPar) {

  }
  void RADAU5DAEIntegrator::plot(int* nr, double* told, double* t, double* z, double* cont, int* lrc, int* n, double* rPar, int* iPar, int* irtrn) {
    double tReal, dt;
    double qSize = system->getqSize();
    double uSize = system->getuSize();
    double H = rPar[0];
//    bool scaling = iPar[4];

    while (*t >= tPlot) {
      for (int i = 1; i <= *n; i++)
        YInp(i - 1) = CONTR5(&i, &tPlot, cont, lrc);

      // H = 1 for unscaling case
      tReal = tPlot * H;
      YInp(qSize, qSize + uSize - 1) = YInp(qSize, qSize + uSize - 1) / H;
      dt = (*t - *told) * H;

      system->plot(YInp, tReal);

      double s1 = clock();
      time += (s1 - s0) / CLOCKS_PER_SEC;
      s0 = s1;

      if (output_)
        cout << "   t = " << tReal << ",\tdt = " << dt << "\r" << flush;

      integPlot << tReal << " " << dt << " " << time << endl;

      tPlot += dtOut;
    }
  }

  void RADAU5DAEIntegrator::integratorSetups(MBSim::DynamicSystemSolver& system_) {
    if (scaling) {
      assert(H > 0);
      tEnd = tEnd / H;
      dtPlot = dtPlot / H;
      dt0 = dt0 / H;
      if(DAEIndex != 3){
        cout << "When scaling technique is applied, the DAEIndex has to be 3. The current DAEIndex = " << DAEIndex << endl;
        throw 1;
      }
    }else {
      H = 1;
      rho = 0.;
    }

    zSize = system_.getzSize();
    laBiSize = system_.getlaSize(); // todo: check whether need to different laBiSize and laUniSize
    YSize = zSize + laBiSize;
    //    if (DAEIndex==21) YSize += laBiSize;
    t = tStart;

    Y.resize(YSize, INIT, 0.0);
    z >> Y(0, zSize - 1);

    if (z0.size())
      z = z0;
    else
      system_.initz(z);

    // aTol and rTol can be chosen both scalars or else both arrays of length NEQ.
    assert(aTol.size() == rTol.size());

    if (aTol.size() == 1)
      iTol = 0; // Scalar
    else {
      iTol = 1; // Vector
      assert(aTol.size() == YSize);
    }

    iout = 1; // Subroutine for output

    rPar.resize(2);
    iPar.resize(5);
    rPar(0) = H;
    rPar(1) = rho;
    iPar(0) = DAEIndex;
    iPar(1) = zSize;
    iPar(2) = YSize;
    iPar(3) = laBiSize;
    iPar(4) = scaling;
    iJac = 0;                      // 0: Jacobi Matrix dF/dY wird über finite Differenzen bestimmt; 1: über Subroutine
    if (useExternalJac)
      iJac = 1;
    mlJac = YSize;                  // N: vollbesetzte Jacobi; <N: Jacobi hat Bandstruktur
    muJac = 0;                      // oberer Bandweite der Jacobi; nicht von Bedeutung falls mlJac = N
    iMas = 1;                      // 0:M =EYE;    1:Subroutine fuer M: MAS(N,AM,LMAS,RPAR,IPAR)
    mlMas = 0;                      // lower Bandwith
    muMas = 0;                      // upper Bandwith

    int lJac = YSize;                   // zur Berechnung der Laenge des Working Space
    int lMas = mlMas + muMas + 1;
    int le = YSize;

    lrWork = (YSize * (lJac + lMas + 3 * le + 12) + 20);
    liWork = 3 * YSize + 20;
    rWork.resize(lrWork, INIT, 0.0);
    iWork.resize(liWork, INIT, 0);

    rWork(0) = 0;                // 0 represent default value 1e-16, rounding unit
    if (dtMax)
      rWork(6) = dtMax;

    /**    DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
     *      INCREASE rWork(2), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS
     *     ARE COSTLY. FOR SMALL SYSTEMS iWork(3) SHOULD BE SMALLER
     *     (0.001D0, SAY). NEGATIV iWork(3) FORCES THE CODE TO
     *     COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP.
     *    DEFAULT 0.001D0.
     */
    rWork(2) = 0.001;  // default for small system.
    if (updateJacobian)
      rWork(2) = -1;

    iWork(0) = 0;                                       // iWork(0) .ne. 0, THE CODE TRANSFORMS THE JACOBIAN MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN.
    iWork(1) = maxSteps;                                  //Maximum Step Numbers
    //iWork(2)                                          //max number of Newton Iterations for the solution of the implizi system; def. 7
    // todo: check!
    iWork(4) = system_.getqSize() + system_.getxSize(); // Dimension of Index 1 Variables (for ODE = N)
    iWork(5) = system_.getuSize();                      // Dimension of Index 2 Variables
    iWork(6) = system_.getlaSize();                                      // Dimension of Index 3 Variables
    //    if (DAEIndex == 21) iWork(6) = 2*laBiSize;
  }

  void RADAU5DAEIntegrator::preIntegrate(MBSim::DynamicSystemSolver& system_) {
    system = &system_;

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
      integPlot << "#3 calculation time [s]:" << endl;
    }
  }

  void RADAU5DAEIntegrator::subIntegrate(MBSim::DynamicSystemSolver& system_, double tEnd) {
    tPlot = t + dtPlot;
    dtOut = dtPlot;
    system_.plot(z, t);

    YInp.resize(YSize);
    output_ = output;

    s0 = clock();

    // call the external Fortran integrator
    RADAU5(&YSize, fdae, &t, Y(), &tEnd, &dt0, rTol(), aTol(), &iTol, jac, &iJac, &mlJac, &muJac, mdae, &iMas, &mlMas, &muMas, plot, &iout, rWork(), &lrWork, iWork(), &liWork, rPar(), iPar(), &idid);
  }

  void RADAU5DAEIntegrator::postIntegrate(MBSim::DynamicSystemSolver& system_) {
    cout.unsetf(ios::scientific);
    if (FlagPlotIntegrator)
      integPlot.close();
    if (FlagPlotIntegrationSum) {
      ofstream integSum((name + ".sum").c_str());
      integSum << "Integration time:  " << time << endl;
      integSum << "Integration steps: " << iWork(16) << "  (total " << iWork(15) << "; " << iWork(17) << " rejected due to error test)" << endl;
      integSum << "Function calls :   " << iWork(13) << "(without numerical evaluation of Jacobian)" << endl;
      integSum << "Jacobian calls :   " << iWork(14) << endl;
      integSum.close();
    }

    if (FlagCoutInfo) {
      if (output)
        cout << endl << endl;
      cout << "Summary Integration with RADAU5DAEIntegrator : " << endl;
      cout << "Integration time:  " << time << endl;
      cout << "Integration steps: " << iWork(16) << "  (total " << iWork(15) << "; " << iWork(17) << " rejected due to error test)" << endl;
      cout << "Function calls :   " << iWork(13) << "(without numerical evaluation of Jacobian)" << endl;
      cout << "Jacobian calls :   " << iWork(14) << endl;
      cout << endl;
    }
  }

}
