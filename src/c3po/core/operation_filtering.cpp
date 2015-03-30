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

#include "operation_filtering.h"
#include "operation_container.h"
#include "selector_container.h"
#include "filter_base.h"
#include "string.h"
#include "comm.h"
#include "output.h"
#include "error.h"
#include <fstream>
#include <cmath>
#include <string>

#include "timer.h"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

OperationFiltering::OperationFiltering(c3po *ptr,const char *name) 
: 
OperationBase(ptr,name),
par_(false)
{

    operationID_ = operationContainer().filterCount(); //set ID, this is the 'old' sampleCount!
}

OperationFiltering::~OperationFiltering()
{
 delete VF_to_filter;
 delete SF_to_filter;
}

// * * * * * * * * * * * * * * * * * * * * * * * *
void OperationFiltering::init(QJsonObject jsonObj) 
{

  if(jsonObj["lagrangian"].isNull())
   error().throw_error_one(FLERR,"You must specify the \"lagrangian\" field!!. \n");
         
  par_=jsonObj["lagrangian"].toBool();
  
  registerInputFields(jsonObj); 
  
  process_input(jsonObj);
  
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
void OperationFiltering::registerInputFields(QJsonObject jsonObj) const
{
   QString qsV=jsonObj["VectorfieldsToFilter"].toString();
   std::string big=qsV.toUtf8().constData();
  
  for (unsigned int it=0; it<big.size();it++)
  {
        
        if (big[it] != ' ' )
        {  
          
          std::string name;
         
          while (big[it] != ' ')
            {
               name.push_back(big[it]);
               it++;  
               if (it==big.size()) break;           
              
            }  
       // if(name.empty()==false)  
        // {
          int index = input().getVFnumber(name);
          if (index==-1) error().throw_error_all("operation_filtering.cpp",0, "Vector fields to filter should be registered in C3PO first");
          filterVFName_.push_back(name);
        // }  
        }
  }
  
  VF_to_filter=new int(filterVFName_.size());
  
  QString qsS=jsonObj["ScalarfieldsToFilter"].toString();
   std::string bigS=qsS.toUtf8().constData();
   
  for (unsigned int it=0; it<bigS.size();it++)
  {
        
        if (bigS[it] != ' ' )
        {  
          
          std::string name;
         
          while (bigS[it] != ' ')
            {
               name.push_back(bigS[it]);
               it++;  
               if (it==bigS.size()) break;           
              
            }  
        //  if(name.empty()==false)  
        // { 
          int index = input().getSFnumber(name);
          if (index==-1) error().throw_error_all("operation_filtering.cpp",0, "Scalar fields to filter should be registered in C3PO first");
          { 
           filterSFName_.push_back(name);
           totalSFName_.push_back(name);
          }
         //} 
        }
  }
  
  SF_to_filter= new int(filterSFName_.size());

}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void OperationFiltering::begin_of_step()
{
 
 for(int i=0;i<*VF_to_filter;i++)
  {
   RMAVF_.push_back(dataStorage().fVFid(filterVFName_[i]));
   if (RMAVF_[i]==-1)
    error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The vector field you want to filter is NOT registered!");
   filterVF_.push_back(dataStorage().VFid(filterVFName_[i]));
    if (filterVF_[i]==-1)
    error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The vector field you want to filter is NOT registered!!");
   std::string name(dataStorage().fVF(filterVF_[i])->name());
    name.append("_");
    name.append(name_);
    dataStorage().fVF(filterVF_[i])->setName(name);
  }
  
  for(int i=0;i<*SF_to_filter;i++)
  {
   RMASF_.push_back(dataStorage().fSFid(filterSFName_[i]));
   if (RMASF_[i]==-1)
    error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The scalar field you want to filter is NOT registered!");
    filterSF_.push_back(dataStorage().SFid(filterSFName_[i]));
    if (filterSF_[i]==-1)
    error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The scalar field you want to filter is NOT registered!!");
    std::string name(dataStorage().fSF(filterSF_[i])->name());
    name.append("_");
    name.append(name_);
    dataStorage().fSF(filterSF_[i])->setName(name);
  }
  
  check_additional_fields();
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void OperationFiltering::setVectorToZero(double* vector, int DIM)
{
 for(int i=0;i<DIM;i++)
  vector[i]=0.0;
}
