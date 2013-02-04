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

#ifndef _CONTACT_FLEXIBLE_BILATERAL_H_
#define _CONTACT_FLEXIBLE_BILATERAL_H_

#include "contact_flexible.h"

namespace MBSim {
  /*! \brief Class for linear flexible bilateral contacts with regularised and decoupled Coulomb friction */
  class ContactFlexibleBilateral: public ContactFlexible {

    public: 
      /*! Constructor */
      ContactFlexibleBilateral(const string &name);
      /*! Destructor */
	  virtual ~ContactFlexibleBilateral() {}
      /*! Evaluates the fully functional force law */
      void updateKinetics(double t);
      /*! Overwrites checkActive of Contact */
      void checkActive() {}
      /*! Compute potential energy */
      double computePotentialEnergy();
  };

}

#endif /* _CONTACT_FLEXIBLE_BILATERAL_H_ */