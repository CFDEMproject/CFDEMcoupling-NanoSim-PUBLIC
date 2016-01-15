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


#include "selector_cellIJK.h"
#include "operation_sampling.h"
#include "operation_container.h"
#include "data_storage.h"
#include "error.h"
#include "mesh.h"
#include "selector_container.h"
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>
#include "input.h"
#include "output.h"
#include "operation_sampling_I_general.h"


using namespace C3PO_NS;


SamplingGeneral::SamplingGeneral(c3po *ptr,const char *name) 
: 
OperationSampling(ptr,name),
component_(-1),
selective_(false)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
SamplingGeneral::~SamplingGeneral()
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

void SamplingGeneral::begin_of_step()
{
  (this->*run)();
}

//-----------------------------------------------------------------------------

void SamplingGeneral::process_input(QJsonObject jsonObj)
{

  
  if(jsonObj["VFieldsToSample"].isNull() && jsonObj["SFieldsToSample"].isNull())
   error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR: You must specify valid 'VFieldsToSample' or 'SFieldsToSample' for general sampling \n");
  if(jsonObj["marker"].isNull())
   error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR: You must specify a valid 'marker' for general sampling \n");
  if(sampleCount_ < 0)
    output().write_screen_one("\nWARNING: Will save all Euler cells for 'general' sampling. This might be expensive. \n");

  if(!jsonObj["component"].isNull())
   component_ = jsonObj["component"].toInt();

 
  std::string sampleVF(jsonObj["VFieldsToSample"].toString().toUtf8().constData());
  std::string sampleSF(jsonObj["SFieldsToSample"].toString().toUtf8().constData());
  std::string markers(jsonObj["marker"].toString().toUtf8().constData());
  
  registerInputFields(sampleVF,sampleSF,markers);
  
  NofMarkers_=markers_.size();
         
  if(VFtoSample_.size()>0 && (component_>2 || component_==-1 ) )
   error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR:Invalid vector component! Valid vector components are:\n 0 (= x )\n 1 (= y )\n 2 (= z ) \n");

  //set formula-based sampling
  if(execFormula_) 
  {
        createSampleVectors(VFtoSample_.size(),SFtoSample_.size(),1);

        //set limits
        if(!jsonObj["formulaLimiterNumerator"].isNull())
        {
            QJsonArray limitsTemp=jsonObj["formulaLimiterNumerator"].toArray();
            formula_->setLimiterNumerator( (limitsTemp.at(0)).toDouble(),  
                                           (limitsTemp.at(1)).toDouble()
                                         );
        }
        if(!jsonObj["formulaLimiterDenominator"].isNull())
        {
            QJsonArray limitsTemp=jsonObj["formulaLimiterDenominator"].toArray();
            formula_->setLimiterDenominator( (limitsTemp.at(0)).toDouble(),  
                                             (limitsTemp.at(1)).toDouble()
                                           );
        }
        if(!jsonObj["formulaLimiterSymmetry"].isNull())
        {
            QJsonArray symmetryTemp=jsonObj["formulaLimiterSymmetry"].toArray();
            formula_->setLimiterSymmetry( (symmetryTemp.at(0)).toBool(), (symmetryTemp.at(1)).toBool() );
            if(   (symmetryTemp.at(0)).toBool() 
               && ( (formula_->getLimitNumerator())[0]<0 || (formula_->getLimitNumerator())[1]<0 ) 
              )
               error().throw_error_one("operation_sampling.cpp",
                                       0,
                                       "ERROR: you cannot set use a formulaLimiterSymmetry[0]=true and have a negative limit for formulaLimiterNumerator! \n");
            if(   (symmetryTemp.at(1)).toBool() 
               && ( (formula_->getLimitDenominator())[0]<0 || (formula_->getLimitDenominator())[1]<0 ) 
              )
               error().throw_error_one("operation_sampling.cpp",
                                       0,
                                       "ERROR: you cannot set use a formulaLimiterSymmetry[1]=true and have a negative limit for formulaLimiterDenominator! \n");
        }

  }
  else createSampleVectors(VFtoSample_.size(),SFtoSample_.size(),0);  

  //Set skips for output
  if(!jsonObj["saveOnlyFormula"].isNull())
    if(jsonObj["saveOnlyFormula"].toBool())
    {
        if(!execFormula_)
           error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR: you cannot set 'saveOnlyFormula' and not provide a formula! \n");
        for(uint iSample=0; iSample<(samples_.size()-1); iSample++)
          sampleSetSkip(iSample);
    }

  run=&SamplingGeneral::sample;
  
  if(!jsonObj["selective"].isNull())
  {
   if(jsonObj["selective"].toBool())
     {
      selective_=true;
      if(jsonObj["max"].isNull())
       error().throw_error_one(FLERR,"You must specify the max domain size when using selective sampling operations. \n"); 
      
      if(jsonObj["min"].isNull())
       error().throw_error_one(FLERR,"You must specify the min domain size when using selective sampling operations. \n"); 
      
      for(int i=0;i<3;i++)
      {
       maximum_[i]=jsonObj["max"].toArray()[i].toDouble();
       minimum_[i]=jsonObj["min"].toArray()[i].toDouble();
      }
      
     }    
  
  }

  
}
/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */
void SamplingGeneral::sample() 
{
   
  int totcells_;
  double sampleValue, marker1;
  int cellID,id;
      //Checks
 
  if(sampleCount_<0)
        totcells_=mesh().NofCells();
  else
        totcells_=sampleCount_;
  std::vector<int> idMarker_;
  for( int mark_=0;mark_<NofMarkers_;mark_++)
  {
    int idMarker1=dataStorage().fSFid(markers_[mark_]);
    if(idMarker1==-1) 
     error().throw_error_one("OperationSampling::sample()",0,"ERROR: Marker1 field not registered!!"); 
   idMarker_.push_back(idMarker1);
  }       
  

 //process Vector Fields
  for(unsigned int samp_=0;samp_<VFtoSample_.size();samp_++)
  {
   id=dataStorage().fVFid(VFtoSample_[samp_]);
   if(id==-1) 
    error().throw_error_one("OperationSampling::sample()",0,"ERROR: Vector field not registered!!"); 
 
 
    for (int cell=0;cell<totcells_;cell++)
    { 
      double markerValue[NofMarkers_];
       cellID=cell;    
              
       if(cellID>=mesh().NofCells()) 
            error().throw_error_one("OperationSampling::sample()",
                                    0,"ERROR: cellID out of range!!"); 
       if(selective_)
       {
        bool skip_=false;
        for(int i=0;i<3;i++)
         if(*(mesh().CellCentre(cell)[i])>maximum_[i] || *(mesh().CellCentre(cell)[i])<minimum_[i] )
          skip_=true;
        
        if(skip_) continue;
       }
         
        sampleValue  = dataStorage().VF(id, component_,cellID);
      
        for( int mark_=0;mark_<NofMarkers_;mark_++)  
           markerValue[mark_]=dataStorage().SF(idMarker_[mark_],cellID);

        insertSample(sampleValue,&markerValue[0],samp_); //pushes the sample into the container
        if(input().verbose())
        {
            std::cout<< "VectorField:sampleValue: "<< sampleValue;
            for(int mark_=0;mark_<NofMarkers_;mark_++)
            {
             std::cout<< "   marker" << mark_ << ": " << markerValue[mark_];
            }
            std::cout<< std::endl;
        }
    }
   }
   
   
 //process Scalar Fields
  for(unsigned int samp_=0;samp_<SFtoSample_.size();samp_++)
  {
   id=dataStorage().fSFid(SFtoSample_[samp_]);
   if(id==-1) 
    error().throw_error_one("OperationSampling::sample()",0,"ERROR: Scalar field not registered!!"); 
   
    for (int cell=0;cell<totcells_;cell++)
    { 
       double markerValue[NofMarkers_];
       cellID=cell;    
       
       if(selective_)
       {
        bool skip_=false;
        for(int i=0;i<3;i++)
         if(*(mesh().CellCentre(cell)[i])>maximum_[i] || *(mesh().CellCentre(cell)[i])<minimum_[i] )
          skip_=true;
        
        if(skip_) continue;
       }
              
       if(cellID>=mesh().NofCells()) 
            error().throw_error_one("OperationSampling::sample()",
                                    0,"ERROR: cellID out of range!!"); 
         
        sampleValue  = dataStorage().SF(id,cellID);
        
        for( int mark_=0;mark_<NofMarkers_;mark_++)  
            markerValue[mark_]=dataStorage().SF(idMarker_[mark_],cellID);

        insertSample(sampleValue,&markerValue[0],VFtoSample_.size()+samp_); //pushes the sample into the container
        if(input().verbose())
        {
            std::cout<< "ScalarField:sampleValue: "<< sampleValue;
            for(int mark_=0;mark_<NofMarkers_;mark_++)
            {
             std::cout<< "   marker" << mark_ << ": " << markerValue[mark_];
            }
            std::cout<< std::endl;
        }
    }
  }  
}    


