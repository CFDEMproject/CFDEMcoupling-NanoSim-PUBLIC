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

#ifdef OPERATION_CLASS

OperationStyle(filtering_FavreRunningVariance, FilteringFavreRunningVariance)

#else


#ifndef C3PO_OPERATION_FILTERING_FAVRERUNNINGVARIANCE_I_H
#define C3PO_OPERATION_FILTERING_FAVRERUNNINGVARIANCE_I_H

#include "operation_filtering.h"
#include <vector>

namespace C3PO_NS
{
 class FilteringFavreRunningVariance : public OperationFiltering
 {
  public:
  
  FilteringFavreRunningVariance(c3po *ptr,const char *_name);
  ~FilteringFavreRunningVariance();

//      void init(int narg, char const* const* arg) {};

  void middle_of_step();
  void end_of_step();
  
  void check_additional_fields();
  
  void process_input(QJsonObject jsonObj);
  
  private:
  
  std::vector<double*>  valueContainerScalars_;
  std::vector<double*>  valueContainerVectors_;
  
  void (FilteringFavreRunningVariance::*run)();
  void (FilteringFavreRunningVariance::*end)();
  
  void runFavreRunningVariance();
  void endFavreRunningVariance();
  void runParFavreRunningVariance();
  void white_end() {}; 
  
  mutable int alpha_;
  mutable int alphaRMA_;  
  mutable std::string phaseFractionFieldName_;
  
  double getAlpha(int selCell);
  double getAlpha_inverted(int selCell);
     
  double (FilteringFavreRunningVariance::*getAlphaValue)(int selCell);
 
  mutable int totFields_; 
  mutable double** fieldsPerProc_;
  mutable double** tmpData_;
  
  mutable double** VTmp_;
  mutable double** Var_;
  mutable int totVar_;
       
 };
}

#endif
#endif
