/* Copyright (C) 2004-2006  Martin Förg
 
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
#define FMATVEC_NO_BOUNDS_CHECK
#include "contact_rigid_bilateral.h"
#include "utils/nonsmooth_algebra.h"
#include "multi_body_system.h"

namespace MBSim {

  ContactRigidBilateral::ContactRigidBilateral(const string &name) : ContactRigid(name) {
    active = true;
  }

  void ContactRigidBilateral::projectGS(double dt) {
    double *a = mbs->getGs()();
    int *ia = mbs->getGs().Ip();
    int *ja = mbs->getGs().Jp();
    gdn(0) = s(0);
    for(int j=ia[laInd]; j<ia[laInd+1]; j++) gdn(0) += a[j]*mbs->getla()(ja[j]);

    la(0) -= rFactor(0)*gdn(0);

    for(int i=1; i<=nFric; i++) {
      gdn(i) = s(i);
      for(int j=ia[laInd+i]; j<ia[laInd+1+i]; j++) gdn(i) += a[j]*mbs->getla()(ja[j]);
    }
    if(nFric==1) la(1) = proxCT2D(la(1)-rFactor(1)*gdn(1),((*mue_fun)(abs(gdn(1))))(0)*fabs(la(0)));
    else if(nFric == 2) la(1,2) = proxCT3D(la(1,2)-rFactor(1)*gdn(1,2),((*mue_fun)(nrm2(gdn(1,2))))(0)*fabs(la(0)));
  }

  void ContactRigidBilateral::solveGS(double dt) {
    assert(nFric <= 1);
    double *a = mbs->getGs()();
    int *ia = mbs->getGs().Ip();
    int *ja = mbs->getGs().Jp();
    gdn(0) = s(0);
    for(int j=ia[laInd]+1; j<ia[laInd+1]; j++) gdn(0) += a[j]*mbs->getla()(ja[j]);

    la(0) = -gdn(0)/a[ia[laInd]];

    if(nFric) {
      gdn(1) = s(1);
      for(int j=ia[laInd+1]+1; j<ia[laInd+2]; j++) gdn(1) += a[j]*mbs->getla()(ja[j]);

      double laNmue = fabs(la(0))*((*mue_fun)(abs(gdn(1))))(0);
      double sdG = -gdn(1)/a[ia[laInd+1]];

      if(fabs(sdG)<=laNmue) la(1) = sdG;
      else la(1) = (laNmue<=sdG) ? laNmue : -laNmue;
    }
  }

  void ContactRigidBilateral::residualProj(double dt) {
    double *a = mbs->getGs()();
    int *ia = mbs->getGs().Ip();
    int *ja = mbs->getGs().Jp();
    for(int i=0; i < 1+nFric; i++) {
      gdn(i) = s(i);
      for(int j=ia[laInd+i]; j<ia[laInd+1+i]; j++) gdn(i) += a[j]*mbs->getla()(ja[j]);
    }

    res(0) = gdn(0);

    if(nFric==1) {
      updateFrictionCoefficient(abs(gdn(1)));
      argT(0) = la(1) - rFactor(1)*gdn(1);
      res(1) = la(1) - proxCT2D(argT(0),mue*fabs(la(0)));

    } else if(nFric == 2) {
      updateFrictionCoefficient(nrm2(gdn(1,2)));
      argT = la(1,2) - rFactor(1)*gdn(1,2);
      res(1,2) = la(1,2) - proxCT3D(argT,mue*fabs(la(0)));
    }
  }

  void ContactRigidBilateral::residualProjJac(double dt) {
    SqrMat Jprox = mbs->getJprox();
    SymMat G = mbs->getG();

    RowVec jp1=Jprox.row(laInd);
    for(int i=0; i<G.size(); i++)
      jp1(i) = G(laInd,i);

    if(nFric == 1) {
      RowVec jp2=Jprox.row(laInd+1);
      double rfac1 = rFactor(1);
      double laNmue = fabs(la(0))*mue;
      if(fabs(argT(0))<=laNmue) for(int i=0; i<G.size(); i++) jp2(i) = rfac1*G(laInd+1,i);
      else {
		jp2.init(0);
		jp2(laInd+1) = 1;
		if(argT(0)>=0) {
		  if(la(0) >= 0) jp2(laInd) = -mue;
		  else jp2(laInd) = mue;
		}
		else {
		  if(la(0) >= 0) jp2(laInd) = mue;
		  else jp2(laInd) = -mue;
		}
      }
    }
    else if(nFric == 2) {
      RowVec jp2=Jprox.row(laInd+1);
      RowVec jp3=Jprox.row(laInd+2);
      double LaT = pow(argT(0),2)+pow(argT(1),2);
      double fabsLaT = sqrt(LaT);
      double laNmue = fabs(la(0))*mue;
      double rFac1 = rFactor(1);
      if(fabsLaT <=  laNmue) {
		for(int i=0; i<G.size(); i++) {
		  jp2(i) = rFac1*G(laInd+1,i);
		  jp3(i) = rFac1*G(laInd+2,i);
		}
      }
      else {
		SymMat dfda(2,NONINIT);
		dfda(0,0) = 1-argT(0)*argT(0)/LaT;
		dfda(1,1) = 1-argT(1)*argT(1)/LaT;
		dfda(0,1) = -argT(0)*argT(1)/LaT;
	
		for(int i=0; i<G.size(); i++) {
		  double e1;
		  if(i==laInd)
		    if(la(0)>=0) e1 = mue;
		    else e1 = -mue;
		  else e1 = 0;
		  double e2 = (i==laInd+1?1.0:0.0);
		  double e3 = (i==laInd+2?1.0:0.0);
		  jp2(i) = e2 - ((dfda(0,0)*(e2 - rFac1*G(laInd+1,i)) + dfda(0,1)*(e3 - rFac1*G(laInd+2,i)))*laNmue + e1*argT(0))/fabsLaT;
		  jp3(i) = e3 - ((dfda(1,0)*(e2 - rFac1*G(laInd+1,i)) + dfda(1,1)*(e3 - rFac1*G(laInd+2,i)))*laNmue + e1*argT(1))/fabsLaT;
		}
      }
    }
  }

  void ContactRigidBilateral::checkForTermination(double dt) {

    double *a = mbs->getGs()();
    int *ia = mbs->getGs().Ip();
    int *ja = mbs->getGs().Jp();
    for(int i=0; i < 1+nFric; i++) {
      gdn(i) = s(i);
      for(int j=ia[laInd+i]; j<ia[laInd+1+i]; j++) gdn(i) += a[j]*mbs->getla()(ja[j]);
    }

    if(fabs(gdn(0)) <= gdTol) ;
    else {
      mbs->setTermination(false);
      return;
    }

    if(nFric==1) {

      if(fabs(gdn(1)) > 0 && fabs(la(1) + gdn(1)/fabs(gdn(1))*((*mue_fun)(abs(gdn(1))))(0)*fabs(la(0))) <= laTol*dt) ;
      else if(fabs(la(1)) <= ((*mue_fun)(abs(gdn(1))))(0)*fabs(la(0)) && fabs(gdn(1)) <= gdTol) ;
      else {
		mbs->setTermination(false);
		return;
      }

    }
    else if(nFric==2) {

      if(nrm2(gdn(1,2)) > 0 && nrm2(la(1,2) + gdn(1,2)/nrm2(gdn(1,2))*((*mue_fun)(nrm2(gdn(1,2))))(0)*fabs(la(0))) <= laTol*dt) ;
      else if(nrm2(la(1,2)) <= ((*mue_fun)(nrm2(gdn(1,2))))(0)*fabs(la(0)) && nrm2(gdn(1,2)) <= gdTol) ;
      else {
		mbs->setTermination(false);
		return;
      }
    }
  }

}