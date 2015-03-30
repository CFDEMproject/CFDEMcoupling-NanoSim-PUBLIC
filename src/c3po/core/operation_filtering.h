/*-----------------------------------------------------------------------------*\
                  ___   _____   _____   _____     ___   
                 / ___\/\  __`\/\  __`\/\  __`\  / __`\ 
                /\ \__/\ \ \_\ \ \ \_\ \ \ \_\ \/\ \_\ \
                \ \____\\ \  __/\ \  __/\ \  __/\ \____/
                 \/____/ \ \ \/  \ \ \/  \ \ \/  \/___/ 
                          \ \_\   \ \_\   \ \_\         
                           \/_/    \/_/    \/_/         

         A Compilation for Fluid-Particle Data Post PrOcessing

Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
               2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria
---------------------------------------------------------------------------------
License
    CPPPO is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with CPPPO. If not, see <http://www.gnu.org/licenses/lgpl.html>.

	This code is designed for on-the-fly post processing of fluid-particle
	data (e.g., of velocity, pressure, concentration, temperature field).

	Parts of the code were developed in the frame of the NanoSim project funded
	by the European Commission through FP7 Grant agreement no. 604656.
\*-----------------------------------------------------------------------------*/

#ifndef C3PO_OPERATION_FILTERING_H
#define C3PO_OPERATION_FILTERING_H

#include "operation_base.h"
#include <vector>

namespace C3PO_NS
{


class OperationFiltering : public OperationBase 
{
  public:

  OperationFiltering(c3po *ptr,const char *_name);
  ~OperationFiltering();

  void begin_of_step();
  virtual void middle_of_step() {};
  virtual void end_of_step() {};
  
  void init(QJsonObject jsonObj);
  void registerInputFields(QJsonObject jsonObj) const;
      
  virtual void process_input(QJsonObject jsonObj) {};
    
  bool particleBased() {return par_;}; 
      
  std::vector<std::string> returnVectorFTFNames() {return filterVFName_;};
  std::vector<std::string> returnScalarFTFNames() {return totalSFName_; };       
      
  protected:
  
  void setVectorToZero(double* vector, int DIM);
  
  virtual void check_additional_fields() {};
    
  mutable std::vector<int> filterVF_;
  mutable std::vector<int> filterSF_;
   
  mutable std::vector<int> RMAVF_;
  mutable std::vector<int> RMASF_;
    
  mutable std::vector<std::string> filterVFName_;
  mutable std::vector<std::string> filterSFName_;
  mutable std::vector<std::string> totalSFName_;
    
  mutable int* VF_to_filter;
  mutable int* SF_to_filter; 
  mutable bool par_;
        
};

} //end c3po_NS

#endif
