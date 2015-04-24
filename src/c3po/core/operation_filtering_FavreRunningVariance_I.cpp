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
#include "operation_filtering_FavreRunningVariance_I.h"
#include "operation_container.h"
#include "selector_container.h"
#include "filter_base.h"
#include "string.h"
#include "mesh.h"
#include "comm.h"
#include "output.h"
#include "error.h"
#include <fstream>
#include <cmath>
#include <string>

#include "timer.h"

#define TOLERANCE_ALPHAVAL 1e-10

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

FilteringFavreRunningVariance::FilteringFavreRunningVariance(c3po *ptr,const char *name) 
: 
OperationFiltering(ptr,name),
alpha_(-1)
{
}

FilteringFavreRunningVariance::~FilteringFavreRunningVariance()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

void FilteringFavreRunningVariance::process_input(QJsonObject jsonObj)
{
  if(jsonObj["phaseFractionField"].isNull())
   error().throw_error_one("operation_filtering.cpp",0,"ERROR: You must edit the \"phaseFractionField\" field in the Json file when using FavreRunningVariance averaging \n");
        
  QString Qalpha=jsonObj["phaseFractionField"].toString();
  std::string alpha=Qalpha.toUtf8().constData();
        
  phaseFractionFieldName_.assign(alpha);
  totalSFName_.push_back(phaseFractionFieldName_);
        
  if(jsonObj["invertPhaseFraction"].toBool())
   getAlphaValue=&FilteringFavreRunningVariance::getAlpha_inverted;
  else
   getAlphaValue=&FilteringFavreRunningVariance::getAlpha;
              
  if(jsonObj["lagrangian"].toBool())
  {
   run=&FilteringFavreRunningVariance::runParFavreRunningVariance;
   end=&FilteringFavreRunningVariance::white_end;
  }
  else
  {
   run=&FilteringFavreRunningVariance::runFavreRunningVariance;
   end=&FilteringFavreRunningVariance::endFavreRunningVariance; 
  }      
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void FilteringFavreRunningVariance::check_additional_fields()
{
  alphaRMA_=dataStorage().fSFid(phaseFractionFieldName_);
  if (alphaRMA_==-1)
   error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The \"phaseFractionField\" must be a scalar field registerd in C3PO \n");
  alpha_=dataStorage().SFid(phaseFractionFieldName_);
  if (alpha_==-1)
   error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The \"phaseFractionField\" must be a scalar field registerd in C3PO. \n");
  std::string name(dataStorage().fSF(alpha_)->name());
  name.append("_");
  name.append(name_);
  dataStorage().fSF(alpha_)->setName(name);
 
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void FilteringFavreRunningVariance::middle_of_step()
{
 (this->*run)();
}

/* * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
void FilteringFavreRunningVariance::end_of_step()
{
 (this->*end)();
}


/* * * * * * * * * * * * * * * * * * * FavreRunningVariance Averaging * * * * * * * * * * * * * * * * * * * * */ 

void FilteringFavreRunningVariance::runFavreRunningVariance()
{
  timer().stamp();

  //define some stuff
  int    selCell=*(selectorContainer().currentCell());
  int    nprocs_=comm().nprocs();
  int*   key=selectorContainer().getKey();

  double alphaValues[nprocs_];
  double alphaValue;
  double value;
  double valueVec[3];

  //tmp fields needed for variance calculations  
  double meanOldStorage[*SF_to_filter];
  double meanNewStorage[*SF_to_filter];
  double meanOldStorageVec[*VF_to_filter][3];
  double meanNewStorageVec[*VF_to_filter][3];

  std::vector<int>* cellList=selectorContainer().getCellsInFilter();
  
  if(selCell >= mesh().MaxNofCellsProc()[comm().me()]) selCell=-1;
  
  //Fill temporary containers with data
  if(selCell != -1) alphaValue = *(mesh().CellVol(selCell))*(this->*getAlphaValue)(selCell);  
  else alphaValue = 0.0;
  MPI_Allgather(&alphaValue, 1, MPI_DOUBLE, &alphaValues, 1, MPI_DOUBLE,MPI_COMM_WORLD);

  valueContainerScalars_.clear();
  for (int i=0;i<*SF_to_filter;i++)
  {   
    if(selCell!=-1)     value = dataStorage().SF(RMASF_[i],selCell);  
    else                value = 0.0;

    valueContainerScalars_.push_back(new double[nprocs_]); //add a new double vector            
    MPI_Allgather(&value, 1, MPI_DOUBLE, valueContainerScalars_[i], 1, MPI_DOUBLE,MPI_COMM_WORLD);
  }
  
  valueContainerVectors_.clear();
  for (int i=0;i<*VF_to_filter;i++)
  {   
    if(selCell!=-1) for (int coord=0;coord<3;coord++)       valueVec[coord] = dataStorage().VF(RMAVF_[i],coord,selCell);
    else            for (int coord=0;coord<3;coord++)       valueVec[coord] = 0.0;

    valueContainerVectors_.push_back(new double[nprocs_*3]); //add a new double vector            
    MPI_Allgather(&valueVec, 3, MPI_DOUBLE, valueContainerVectors_[i], 3, MPI_DOUBLE,MPI_COMM_WORLD);
  }
  
  timer().stamp(TIME_FILTER_MPI);
  timer().stamp();

  //Perform main filtering loop
  int point=0;
  for (int k=0;k<nprocs_;k++)
  { 
        for (int j=point;j<key[k]+point;j++) 
        {
            alphaValue = alphaValues[k]; 
            dataStorage().fSF(alpha_)->value()[(*cellList)[j]] += alphaValue;
            double Wn = dataStorage().fSF(alpha_)->value()[(*cellList)[j]] ;
            double alphaDivW = alphaValue
                                / (
                                      Wn 
                                    + *(mesh().CellVol(selCell)) * TOLERANCE_ALPHAVAL 
                                  );

                                      
            //Loop scalar fields to filter and do variance calculation
            for (int i=0;i<*SF_to_filter;i++)
            { 
                 meanOldStorage[i] = dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]];
                 meanNewStorage[i] = meanOldStorage[i] 
                                   + alphaDivW 
                                   * (valueContainerScalars_[i][k] - meanOldStorage[i]);
                         
                dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]] = meanNewStorage[i];
            }

            //Loop vector fields
            for (int i=0;i<*VF_to_filter;i++)
            { 
               for (int coord=0;coord<3;coord++)
               {
                 meanOldStorageVec[i][coord] = *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]);
                 meanNewStorageVec[i][coord] = meanOldStorageVec[i][coord]
                         + alphaDivW 
                         * (valueContainerVectors_[i][k*3+coord] - meanOldStorageVec[i][coord]);
                         
                 *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]) = meanNewStorageVec[i][coord];
               }
            }
            
            //Loop Variance fields:
            if(!computeVariance_)            
                continue;


            double varianceOld;
            double varianceNew;

            //Loop SCALAR fields and do variance calculation
            for(int kVar=0;kVar<*SF_for_varianceCalc_;kVar++)
            {
                 int i = *SF_to_filter + kVar;
                 varianceOld = dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]]; //no separate storage, uses filter storage!

                 if(varianceHasCrossTerm_)
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                    valueContainerScalars_[filterSFVarianceValueID_[kVar]][k]       //xi,n
                                   *valueContainerScalars_[filterSFVarianceValueSecondID_[kVar]][k] //yi,n
                                  - meanOldStorage[filterSFVarianceValueID_[kVar]]                  //\overbar{xi,n-1]
                                   *meanOldStorage[filterSFVarianceValueSecondID_[kVar]]            //\overbar{yi,n-1]
                                 )
                               + Wn
                               *(
                                    meanOldStorage[filterSFVarianceValueID_[kVar]]                  //\overbar{xi,n-1]
                                   *meanOldStorage[filterSFVarianceValueSecondID_[kVar]]            //\overbar{yi,n-1]
                                  - meanNewStorage[filterSFVarianceValueID_[kVar]]                  //\overbar{xi,n]
                                   *meanNewStorage[filterSFVarianceValueSecondID_[kVar]]            //\overbar{yi,n]
                                );
                 else
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                  valueContainerScalars_[filterSFVarianceValueID_[kVar]][k] 
                                - meanOldStorage[filterSFVarianceValueID_[kVar]]
                                 )
                               * (
                                  valueContainerScalars_[filterSFVarianceValueID_[kVar]][k] 
                                - meanNewStorage[filterSFVarianceValueID_[kVar]]
                                 );
                         
                dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]] = varianceNew;
            }


            //Loop VECTOR fields and do variance calculation
            for(int kVar=0;kVar<*VF_for_varianceCalc_;kVar++)
            {
               int firstCoord[3]  = {0,1,2};
               int secondCoord[3] = {0,1,2};
               if(filterVFVarianceComputeOffDiagonal_[kVar])
               {
                    firstCoord [1]=0;firstCoord [2]=1;
                    secondCoord[0]=1;secondCoord[1]=2;
               }
               int i = *VF_to_filter + kVar;
               for (int coord=0;coord<3;coord++)
               {

                 varianceOld  = *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]);
    
                 if(evaluateVarianceVectorScalarMixed_[kVar])
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                    valueContainerVectors_[filterVFVarianceValueID_[kVar]][k*3+coord]                 //xi,n
                                   *valueContainerScalars_[filterSFVarianceValueVectorScalarMixed_[kVar]][k]          //yi,n
                                  - meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]                          //\overbar{xi,n-1]
                                   *meanOldStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]                    //\overbar{yi,n-1]
                                 )
                               + Wn
                               * (
                                    meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]                    //\overbar{xi,n-1]
                                   *meanOldStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]             //\overbar{yi,n-1]
                                  - meanNewStorageVec[filterVFVarianceValueID_[kVar]][coord]                   //\overbar{xi,n]
                                   *meanNewStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]            //\overbar{yi,n]
                                 );

                 else if(varianceHasCrossTerm_ || filterVFVarianceComputeOffDiagonal_[kVar] )
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                    valueContainerVectors_[filterVFVarianceValueID_[kVar]][k*3+firstCoord[coord]]          //xi,n
                                   *valueContainerVectors_[filterVFVarianceValueSecondID_[kVar]][k*3+secondCoord[coord]]   //yi,n
                                  - meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]                   //\overbar{xi,n-1]
                                   *meanOldStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]            //\overbar{yi,n-1]
                                 )
                               + Wn
                               * (
                                    meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]                    //\overbar{xi,n-1]
                                   *meanOldStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]             //\overbar{yi,n-1]
                                  - meanNewStorageVec[filterVFVarianceValueID_[kVar]][coord]                   //\overbar{xi,n]
                                   *meanNewStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]            //\overbar{yi,n]
                                 );
                 else
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                  valueContainerVectors_[filterVFVarianceValueID_[kVar]][k*3+coord]
                                - meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]
                                 )
                               * (
                                  valueContainerVectors_[filterVFVarianceValueID_[kVar]][k*3+coord] 
                                - meanNewStorageVec[filterVFVarianceValueID_[kVar]][coord]
                                 );
                         
                 *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]) = varianceNew;

//                  if(coord==0)
//                  {
//                     std::cout << "varianceNew: " << varianceNew << ", varianceOld: " << varianceOld 
//                               << ", alphaValue: " << alphaValue << ", Wn: " << Wn 
//                               << ", value1: " << valueContainerVectors_[filterVFVarianceValueID_[kVar]][k*3+firstCoord[coord]] 
//                               << ", value2: " << valueContainerVectors_[filterVFVarianceValueSecondID_[kVar]][k*3+secondCoord[coord]]  
//                               << ", meanOldStorageVec[" << filterVFVarianceValueID_[kVar] << "]: " << meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]
//                               << ", meanNewStorageVec[" << filterVFVarianceValueID_[kVar] << "]: " << meanNewStorageVec[filterVFVarianceValueID_[kVar]][coord]  
//                               << " \n"; 
//                  }
               }
            }
            
        }
        point+=key[k];

  }

  timer().stamp(TIME_FILTER);
// END MAIN FILTER LOOP
          
}
/*---------------------------------------------------------------------*/
void FilteringFavreRunningVariance::endFavreRunningVariance()
{

  for (int cell=0;cell<mesh().NofCells();cell++) 
  {
     double alphaValue=dataStorage().fSF(alpha_)->value()[cell];
     dataStorage().fSF(alpha_)->value()[cell] /= (*( selectorContainer().filterVolume(cell)) );

     if(!computeVariance_)
        continue;

     //Finalize variance calculation - SCALARS
     for(int kVar=0;kVar<*SF_for_varianceCalc_;kVar++)
     {
          int i = *SF_to_filter + kVar;
          if(abs(alphaValue)>(TOLERANCE_ALPHAVAL*(*( selectorContainer().filterVolume(cell)) ))) //must multply with cell volume!
            for (int coord=0;coord<3;coord++)
               dataStorage().fSF(filterSF_[i])->value()[cell] /= alphaValue;
          else
            for (int coord=0;coord<3;coord++)
               dataStorage().fSF(filterSF_[i])->value()[cell] = 0.0;
     }

     //Finalize variance calculation - VECTORS
     for(int kVar=0;kVar<*VF_for_varianceCalc_;kVar++)
     {
          int i = *VF_to_filter + kVar;
          if(abs(alphaValue)>(TOLERANCE_ALPHAVAL*(*( selectorContainer().filterVolume(cell)) ))) //must multply with cell volume!
            for (int coord=0;coord<3;coord++)
               *dataStorage().fVF(filterVF_[i])->value(coord,cell) /= alphaValue;
          else
            for (int coord=0;coord<3;coord++)
               *dataStorage().fVF(filterVF_[i])->value(coord,cell) = 0.0;
     }
  }

}
/*------------------------------------------------------------------------*/
void FilteringFavreRunningVariance::runParFavreRunningVariance()
{
   //TODO: needs to be implemented 
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double FilteringFavreRunningVariance::getAlpha(int selCell)
{ 
 return dataStorage().SF(alphaRMA_,selCell); 
}

double FilteringFavreRunningVariance::getAlpha_inverted(int selCell)
{ 
 return ( 1.0 - dataStorage().SF(alphaRMA_,selCell) );  
}

