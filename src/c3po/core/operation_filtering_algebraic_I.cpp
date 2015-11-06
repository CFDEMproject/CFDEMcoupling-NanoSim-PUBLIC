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
#include "memory_ns.h"

#define TOLERANCE_VOL 1e-10

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
 delete displs_;
 delete recvCount_;
 delete recvBuf_;
 delete fieldsToAdd_;
 
 C3PO_MEMORY_NS::destroy(fieldsPerProc_);
 C3PO_MEMORY_NS::destroy(tmpData_);
 C3PO_MEMORY_NS::destroy(VTmp_);
 C3PO_MEMORY_NS::destroy(Var_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

void FilteringAlgebraic::process_input(QJsonObject jsonObj)
{
  

  if(jsonObj["lagrangian"].toBool())
  {
   run=&FilteringAlgebraic::runParAlgebraic;
   end=&FilteringAlgebraic::white_end;
   
   totFields_=*SF_to_filter + 3*(*VF_to_filter)+1;
   
   fieldsPerProc_= C3PO_MEMORY_NS::create(fieldsPerProc_, comm().nprocs(), totFields_);
   
   
   tmpData_=C3PO_MEMORY_NS::create(tmpData_, comm().nprocs(), totFields_);
   
   if(computeVariance_)            
   {
    totVar_=*SF_for_varianceCalc_+(*VF_for_varianceCalc_*3);
    VTmp_= C3PO_MEMORY_NS::create(VTmp_, comm().nprocs(), totVar_);
    Var_= C3PO_MEMORY_NS::create(Var_, comm().nprocs(), totVar_);
   
   } 
   
  }
  else
  {
   run=&FilteringAlgebraic::runAlgebraic;
   end=&FilteringAlgebraic::endAlgebraic;  
   
   totFields_=*SF_to_filter + 3*(*VF_to_filter);
   
  fieldsPerProc_= C3PO_MEMORY_NS::create(fieldsPerProc_, comm().nprocs(), totFields_);
   
     //prepare for MPI
   displs_=new int[comm().nprocs()];
   recvCount_=new int[comm().nprocs()];
   double buf=0;
   for(int n=0;n<comm().nprocs();n++)
   {
    displs_[n]=buf;
    buf+=totFields_;
    recvCount_[n] = totFields_;
   }
  
  
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
  int* key=selectorContainer().getKey();
  std::vector<int>* cellList=selectorContainer().getCellsInFilter();
  double commFields_[totFields_];
  

 if(selCell < mesh().MaxNofCellsProc()[comm().me()])
 {
  int n=0;
  double vol =*(mesh().CellVol(selCell));
  for (int i=0;i<*VF_to_filter;i++) 
   for (int coord=0;coord<3;coord++)
   {
    commFields_[n]= dataStorage().VF(RMAVF_[i],coord,selCell)*vol;
    n++;
   } 
   
  for (int i=0;i<*SF_to_filter;i++)
  {
   commFields_[n]=dataStorage().SF(RMASF_[i],selCell)*vol; 
   n++;
  } 

 } 
 else
 {
  for(int n=0;n<totFields_;n++)
   commFields_[n]=0;
 }
  
 timer().stamp(TIME_FILTER);
 timer().stamp();
  
 MPI_Allgatherv(&commFields_,totFields_, MPI_DOUBLE, &fieldsPerProc_[0][0], recvCount_, displs_,MPI_DOUBLE, MPI_COMM_WORLD);
  
 timer().stamp(TIME_FILTER_MPI);
 timer().stamp();


//Sum Vector Fields 
int n=0;
for (int i=0;i<*VF_to_filter;i++)
{ 
//    printf("runFavre: proocessing vector field[%d] with id: %d \n", i,filterVF_[i]);

    for (int coord=0;coord<3;coord++)
    {
       
        int point=0;
        for (int k=0;k<nprocs_;k++)
        { 
          
          for (int j=point;j<key[k]+point;j++) 
           {
             
             *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]) += fieldsPerProc_[k][n];
                       
//             printf("coord: %d, process: %d, point: %d, value: %g, pulled value: %g, pulled alpha: %g  \n",
 //                    coord, k, j, *dataStorage().fVF(filterVF_[i])->value(coord,cellList[j]), value,alphaValue );
          }
        point+=key[k];

       }
       n++;    
    } 
}
//Sum Scalar Fields 
for (int i=0;i<*SF_to_filter;i++)
{ 
//    printf("runFavre: proocessing scalar field[%d] with id: %d \n", i,filterSF_[i]);

    int point=0;
   
    for (int k=0;k<nprocs_;k++)
     { 

       for (int j=point;j<key[k]+point;j++) 
       {
         
         dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]] +=fieldsPerProc_[k][n];
//         printf("process: %d, point: %d, value: %g. \n",
//                 k, j, dataStorage().fSF(filterSF_[i])->value()[cellList[j]] );
       }
       point+=key[k];
     }
     n++;    
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
  int selPar=selectorContainer().currentParticle();
  int nprocs_=comm().nprocs();
  std::vector<int>* cellList=selectorContainer().getCellsInFilter();
  int* key = selectorContainer().getKey();
  int totvec= 3*(*VF_to_filter);
  double localAlpha_;
  int point=0;
  
  if(selPar >=  dataStorage().getNofParticlesProc(comm().me())) selPar=-1;
  
  setVectorToZero(&fieldsPerProc_[0][0],totFields_*nprocs_);

  //Calculating partial sum
  for (int k=0;k<nprocs_;k++)
  { 
   
      
    for (int s=point;s<key[k]+point;s++) 
    {
     
     localAlpha_= *(mesh().CellVol((*cellList)[s]));    
     for(int i=0;i<*VF_to_filter;i++)
      for (int coord=0;coord<3;coord++)
       fieldsPerProc_[k][3*i+coord] += localAlpha_ * dataStorage().VF(RMAVF_[i],coord,(*cellList)[s]);
      
    
     for(int i=0;i<*SF_to_filter;i++)  
      fieldsPerProc_[k][i+totvec] += localAlpha_ * dataStorage().SF(RMASF_[i],(*cellList)[s]);
     
     fieldsPerProc_[k][totFields_-1]+= localAlpha_;
     
    }
    

   point+=key[k];

  }
  
  timer().stamp(TIME_FILTER);
  timer().stamp();
  
   MPI_Allreduce(&fieldsPerProc_[0][0],&tmpData_[0][0], totFields_*nprocs_, MPI_DOUBLE, MPI_SUM ,MPI_COMM_WORLD); 
   
  timer().stamp(TIME_FILTER_MPI);
  timer().stamp();
  
 if(selPar>-1 && tmpData_[comm().me()][totFields_-1]>TOLERANCE_VOL) 
 {
 //Finish Calculation
 
   for(int i=0;i<*VF_to_filter;i++)
      for (int coord=0;coord<3;coord++)
        dataStorage().getParticle(selPar)->filteredVector(filterVF_[i])[coord] = tmpData_[comm().me()][3*i+coord]/tmpData_[comm().me()][totFields_-1];
        
  
   for(int i=0;i<*SF_to_filter;i++)  
    *(dataStorage().getParticle(selPar)->filteredScalar(filterSF_[i])) =  tmpData_[comm().me()][i+totvec]/tmpData_[comm().me()][totFields_-1]; 
  }  
       
 timer().stamp(TIME_FILTER);
    
  

    
 
  if(computeVariance_)            
  { 
    point=0;
    
    setVectorToZero(&VTmp_[0][0],totVar_*nprocs_);
    int varVecs=(*VF_for_varianceCalc_*3);
    for (int k=0;k<nprocs_;k++)
    {
     if(tmpData_[k][totFields_-1]>TOLERANCE_VOL)
     {
      for (int s=point;s<key[k]+point;s++) 
      {
       localAlpha_= *(mesh().CellVol((*cellList)[s]));
       //Loop SCALAR fields and do variance calculation
        for(int kVar=0;kVar<*SF_for_varianceCalc_;kVar++)
        {
     
         if(varianceHasCrossTerm_)
          VTmp_[k][kVar+varVecs] += localAlpha_                              
                          *( dataStorage().SF(filterSFVarianceValueID_[kVar],(*cellList)[s]) - tmpData_[k][filterSFVarianceValueID_[kVar]+totvec]/tmpData_[k][totFields_-1])
                          *( dataStorage().SF(filterSFVarianceValueSecondID_[kVar],(*cellList)[s]) -  tmpData_[k][filterSFVarianceValueSecondID_[kVar]+totvec]/tmpData_[k][totFields_-1]);
         else
          VTmp_[k][kVar+varVecs] += localAlpha_                             
                          *( dataStorage().SF(filterSFVarianceValueID_[kVar],(*cellList)[s]) - tmpData_[k][filterSFVarianceValueID_[kVar]+totvec]/tmpData_[k][totFields_-1])
                          *( dataStorage().SF(filterSFVarianceValueID_[kVar],(*cellList)[s]) - tmpData_[k][filterSFVarianceValueID_[kVar]+totvec]/tmpData_[k][totFields_-1]);
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
  
         for (int coord=0;coord<3;coord++)
         {

          if(evaluateVarianceVectorScalarMixed_[kVar])
             VTmp_[k][3*kVar+coord] += localAlpha_                              
                                      *( dataStorage().VF(filterVFVarianceValueID_[kVar],coord,(*cellList)[s]) - tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1])
                                     *( dataStorage().VF(filterVFVarianceValueSecondID_[kVar],coord,(*cellList)[s]) - tmpData_[k][3*filterVFVarianceValueSecondID_[kVar]+coord]/tmpData_[k][totFields_-1]);
          else if(varianceHasCrossTerm_ || filterVFVarianceComputeOffDiagonal_[kVar] )
             VTmp_[k][3*kVar+coord] += localAlpha_                              
                                      *( dataStorage().VF(filterVFVarianceValueID_[kVar],firstCoord[coord],(*cellList)[s])  -  tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1])
                                      *( dataStorage().VF(filterVFVarianceValueSecondID_[kVar],secondCoord[coord],(*cellList)[s]) -  tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1]);
          else
             VTmp_[k][3*kVar+coord] += localAlpha_                              
                                      *( dataStorage().VF(filterVFVarianceValueID_[kVar],coord,(*cellList)[s]) -  tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1])
                                      *( dataStorage().VF(filterVFVarianceValueID_[kVar],coord,(*cellList)[s]) -  tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1]);



         }
        }
      }
     }
    }
  
  timer().stamp(TIME_FILTER);
  timer().stamp();
  
   MPI_Allreduce(&VTmp_[0][0],&Var_[0][0], totVar_*nprocs_, MPI_DOUBLE, MPI_SUM ,MPI_COMM_WORLD); 
   
  timer().stamp(TIME_FILTER_MPI);
  timer().stamp();
  if(selPar==-1) return;
  if(tmpData_[comm().me()][totFields_-1]<TOLERANCE_VOL) return;
  
  double finalAlpha_=tmpData_[comm().me()][totFields_-1];
 
   for(int kVar=0;kVar<*VF_for_varianceCalc_;kVar++)
   { 
      int i = *VF_to_filter + kVar;
      for(int coord=0;coord<3;coord++)
      {
        dataStorage().getParticle(selPar)->filteredVector(filterVF_[i])[coord] =  Var_[comm().me()][kVar+coord] /finalAlpha_;     
      }
   }
   
 /*  for(int kVar=0;kVar<*SF_for_varianceCalc_;kVar++)
   { 
      int i = *SF_to_filter + kVar;
       *(dataStorage().getParticle(selPar)->filteredScalar(i)) =  Var_[comm().me()][(*VF_for_varianceCalc_*3)+kVar] /  finalAlpha_; 
   }
  
*/
  
  }
 
   
 
    timer().stamp(TIME_FILTER);
}   
