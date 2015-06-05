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

#ifndef _PHEM56_INTEGRATOR_H_
#define _PHEM56_INTEGRATOR_H_

#include "integratorDAE.h"
#include <fstream>
namespace MBSimIntegrator {

  /** \brief DAE-Integrator PHEM56
   */
  class PHEM56Integrator : public IntegratorDAE {

    private:

      static void fprob(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
      static void plot(int* NSTEP, int* nq, int* nu, int* nx, int* nla, int* LRDO, double* q, double* u, double* x, double* ud, double* RLAM, double* DOWK);
      static void fprob1(int* IFCN, int* qSize, int* uSize, int* xSize, int* laSize, int* LDG, int* LDF, int* LDA, int* NBLK, int* NMRC, int* NPGP, int* NPFL, int* INDGR, int* INDGC, int* INDFLR, int* INDFLC, double* t, double* q_, double* u_, double* x_, double* la_, double* g, double* GQ, double* h, double* GQQ, double* w_hat, double* FL, double* qd, double* xd, double* M);

      static double tPlot;
      static double dtOut;
      static fmatvec::Vec zInp;
      static std::ofstream integPlot;
      static double s0;
      static double time;
      static bool output_;  // redeclare a output flag as it has to be used in a static functio.n
      static Vec rPar;
      static Vector<Var, int> iPar;

      /**
       * \brief linear algebra routines according to the problem formulation:
       *  0: full, L  = G^T, DEC/SOL: for little dimensional system.
       *  1: full, L != G^T, DEC/SOL
       *  2: full, L  = G^T, DGTREF/DGETRS (default): form LAPACK, for general square matrix A and larger dimensional system
       *  3: full, L != G^T, DGTREF/DGETRS
       *  4: sparse, L  = G^T, MA28PACK
       *  5: sparse, L != G^T, MA28PACK
       */
      int mode;

      /**
       * \brief  Additional system state vectors
       */
      fmatvec::Vec q, u, x, ud, la;

      static std::vector<fmatvec::Vec2> la_plot;
      static std::ofstream la0Sum, la1Sum;

    public:

      PHEM56Integrator();
      ~PHEM56Integrator() {}

      void integratorSetups(MBSim::DynamicSystemSolver& system);
      void preIntegrate(MBSim::DynamicSystemSolver& system);
      void subIntegrate(MBSim::DynamicSystemSolver& system, double tEnd);
      void postIntegrate(MBSim::DynamicSystemSolver& system);
  };

}

#endif
