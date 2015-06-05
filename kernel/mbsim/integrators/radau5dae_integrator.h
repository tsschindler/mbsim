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

#ifndef _RADAU5DAE_INTEGRATOR_H_
#define _RADAU5DAE_INTEGRATOR_H_

using namespace fmatvec;

#include "integratorDAE.h"

namespace MBSimIntegrator {

  /** \brief DAE-Integrator RADAU5 (up to index 3)
   */
  class RADAU5DAEIntegrator : public IntegratorDAE {

    private:
      static void fdae(int* YSize, double* t, double* Y_, double* F_, double* rPar, int* iPar);
      static void mdae(int* YSize, double* MasMat_, int* lMas, double* rPar, int* iPar);
      static void plot(int* nr, double* told, double* t, double* z, double* cont, int* lrc, int* n, double* rPar, int* iPar, int* irtrn);
      static void jac(int* YSize, double* t, double* Y_, double* jac_, int* ljac_, double* rpar, int* ipar);

      static double tPlot;
      static double dtOut;
      static fmatvec::Vec YInp;
      static std::ofstream integPlot;
      static double s0;
      static double time;
      static bool output_;

      /**
       * \brief flag of whether update Jacobian matrix at every accepted step.
       */
      bool updateJacobian;

      /**
       * \brief Switch for the computation of the Jacobian. If IJAC=0 the Jacobian is computed internally
       *        by finite differences; the subroutine JAC is not necessary. If IJAC=1 the Jacobian is supplied
       *        by the subroutine JAC.
       */
      int iJac;

      /**
       * \brief Switch for the banded structure of the Jacobian.
       * If mlJac = N: full Sized  Jacobian Matrix;  mlJac < N: Jacobian hat band structure
       */
      int mlJac;

      /**
       * \brief The upper bandwith of the Jacobian matrix.
       * Does not need to be set if mlJac = N;
       */
      int muJac;

      /**
       * \brief Flag of mass matrix. If iMas = 0 then M is take as the identitz matrix.
       * If iMas = 1 the mass matrix is supplied by the user through the subroutine MAS.
       * If iMas = 2 the problem is explicit and neutral.
       */
      int iMas;

      /**
       * \brief Lower bandwidth of the mass matrix
       */
      int mlMas;

      /**
       * \brief Upper bandwidth of the mass matrix
       */
      int muMas;

    public:

      RADAU5DAEIntegrator();
      ~RADAU5DAEIntegrator() {}

      void setDAEIndex(int index);
      void updateJacobianEveryStep(bool flag=true) {updateJacobian = flag;}
//      void setGGLStabilisation() {DAEIndex=21;};

      /** subroutines used by integrate */
      void integratorSetups(MBSim::DynamicSystemSolver& system_);
      void preIntegrate(MBSim::DynamicSystemSolver& system_);
      void subIntegrate(MBSim::DynamicSystemSolver& system_, double tEnd);
      void postIntegrate(MBSim::DynamicSystemSolver& system_);
  };

}

#endif
