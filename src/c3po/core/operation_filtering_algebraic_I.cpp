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
#include "operation_filtering_algebraic_I.h"
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

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

FilteringAlgebraic::FilteringAlgebraic(c3po *ptr,const char *name) 
: 
OperationFiltering(ptr,name) 
{
}

FilteringAlgebraic::~FilteringAlgebraic()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

void FilteringAlgebraic::process_input(QJsonObject jsonObj)
{
  if(jsonObj["lagrangian"].toBool())
  {
   run=&FilteringAlgebraic::runParAlgebraic;
   end=&FilteringAlgebraic::white_end;
  }
  else
  {
   run=&FilteringAlgebraic::runAlgebraic;
   end=&FilteringAlgebraic::endAlgebraic;  
  }
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void FilteringAlgebraic::middle_of_step()
{
 (this->*run)();
}

/* * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
void FilteringAlgebraic::end_of_step()
{
 (this->*end)();
}

/* * * * * * * * * * * * * * * * * * * Algebraic Averaging * * * * * * * * * * * * * * * * * * */
void FilteringAlgebraic::runAlgebraic()
{
 timer().stamp();
  //define some stuff
  int selCell=*(selectorContainer().currentCell());
  int nprocs_=comm().nprocs();
  double values[nprocs_];
  int* key=selectorContainer().getKey();
  double value;
  double alphaValue;
  double alphaValues[nprocs_];
  std::vector<int>* cellList=selectorContainer().getCellsInFilter();
  
  if(selCell >= mesh().MaxNofCellsProc()[comm().me()]) selCell=-1;
  
alphaValue = *(mesh().CellVol(selCell));  
MPI_Allgather(&alphaValue, 1, MPI_DOUBLE, &alphaValues, 1, MPI_DOUBLE,MPI_COMM_WORLD);
  
timer().stamp(TIME_FILTER_MPI);
timer().stamp();

//Sum Vector Fields 
for (int i=0;i<*VF_to_filter;i++)
{ 
//    printf("runFavre: proocessing vector field[%d] with id: %d \n", i,filterVF_[i]);

    for (int coord=0;coord<3;coord++)
    {
       
        if(selCell!=-1) value = dataStorage().VF(RMAVF_[i],coord,selCell);
        else value = 0;
          
        MPI_Allgather(&value, 1, MPI_DOUBLE, &values, 1, MPI_DOUBLE,MPI_COMM_WORLD);
        int point=0;
        for (int k=0;k<nprocs_;k++)
        { 
          for (int j=point;j<key[k]+point;j++) 
           {
             value=values[k];
             alphaValue=alphaValues[k];
             *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]) += value*alphaValue;
             
//             printf("coord: %d, process: %d, point: %d, value: %g, pulled value: %g, pulled alpha: %g  \n",
 //                    coord, k, j, *dataStorage().fVF(filterVF_[i])->value(coord,cellList[j]), value,alphaValue );
        }
        point+=key[k];

     }    
    } 
}
//Sum Scalar Fields 
for (int i=0;i<*SF_to_filter;i++)
{ 
//    printf("runFavre: proocessing scalar field[%d] with id: %d \n", i,filterSF_[i]);

    int point=0;
    if(selCell!=-1) value = dataStorage().SF(RMASF_[i],selCell);  
    else value = 0;
    
    MPI_Allgather(&value, 1, MPI_DOUBLE, &values, 1, MPI_DOUBLE,MPI_COMM_WORLD);
   
    for (int k=0;k<nprocs_;k++)
     { 
       
       for (int j=point;j<key[k]+point;j++) 
       {
         value=values[k];
         alphaValue=alphaValues[k];
         dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]] += value*alphaValue;
//         printf("process: %d, point: %d, value: %g. \n",
//                 k, j, dataStorage().fSF(filterSF_[i])->value()[cellList[j]] );
       }
       point+=key[k];
     }    
} 



timer().stamp(TIME_FILTER);
   //   std::cout << "\nLet's see during filtering..." << *(dataStorage().fVF(0)->value(0,(*cellList)[0])) << " " << (*cellList)[0];
//Calculate FilterSize

//std::cout <<"\n FilterSize: " << fsize_[comm().me()];
           
} 
/*-----------------------------------------------------------------------*/

void FilteringAlgebraic::endAlgebraic()
{
  int NofCells = mesh().NofCells();
  int field_;
  
   for (int cell=0;cell<NofCells;cell++) 
    {
     for (int i=0;i<*VF_to_filter;i++)
     {
      field_=filterVF_[i];
      for (int coord=0;coord<3;coord++)
       *dataStorage().fVF(field_)->value(coord,cell) = *dataStorage().fVF(field_)->value(coord,cell)/(*( selectorContainer().filterVolume(cell)));
     }
     for (int i=0;i<*SF_to_filter;i++)
     {
      field_=filterSF_[i];
      dataStorage().fSF(field_)->value()[cell] =dataStorage().fSF(field_)->value()[cell]/(*( selectorContainer().filterVolume(cell)));
     }
    }

}

/*---------------------------------------------------------------------*/
void FilteringAlgebraic::runParAlgebraic()
{
  timer().stamp();
  //define some stuff
  int selCell=*(selectorContainer().currentCell());
  int nprocs_=comm().nprocs();
  std::vector<int>* cellList=selectorContainer().getCellsInFilter();
  int* key = selectorContainer().getKey();
  double VtmpData_[*VF_to_filter*3]; 
  double StmpData_[*SF_to_filter] ; 
  double StmpAlpha_=0;
 
 
  if(selCell >= mesh().MaxNofCellsProc()[comm().me()]) selCell=-1;
  
  double VData[*VF_to_filter*3];
  double SData[*SF_to_filter ];

  double localVol_;  
 
  int point=0;

  //Calculating partial sum
  for (int k=0;k<nprocs_;k++)
  { 
   setVectorToZero(&VtmpData_[0],*VF_to_filter*3);
   setVectorToZero(&StmpData_[0],*SF_to_filter);
   StmpAlpha_=0.0;
    
    for (int s=point;s<key[k]+point;s++) 
    {
     localVol_= *(mesh().CellVol((*cellList)[s]));
          
     for(int i=0;i<*VF_to_filter;i++)
      for (int coord=0;coord<3;coord++)
       VtmpData_[3*i+coord] += localVol_ * dataStorage().VF(RMAVF_[i],coord,(*cellList)[s]);
      
     for(int i=0;i<*SF_to_filter;i++)  
      StmpData_[i] += localVol_ * dataStorage().SF(RMASF_[i],(*cellList)[s]);
     
     StmpAlpha_ += localVol_;
     
    }
    
  
   MPI_Reduce(VtmpData_,VData, *VF_to_filter*3, MPI_DOUBLE, MPI_SUM ,k,MPI_COMM_WORLD); 
   MPI_Reduce(StmpData_,SData,*SF_to_filter,MPI_DOUBLE,MPI_SUM,k,MPI_COMM_WORLD); 
   
   point+=key[k]; 
  }
  
  
 if(selCell==-1) return;
 
 //Finish Calculation
 
  
   for(int i=0;i<*VF_to_filter;i++)
      for (int coord=0;coord<3;coord++)
        *dataStorage().fVF(filterVF_[i])->value(coord,selCell) = VData[3*i+coord] /  (*( selectorContainer().filterVolume(selectorContainer().filterVolumeSize()-1)) );
        
  
   for(int i=0;i<*SF_to_filter;i++)  
    dataStorage().fSF(filterSF_[i])->value()[selCell] =  SData[i] / (*( selectorContainer().filterVolume(selectorContainer().filterVolumeSize()-1)) ); 
   
       
 
    timer().stamp(TIME_FILTER);
}   
