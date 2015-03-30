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

#include "running_stat.h"
#include <cmath>
#include <sys/stat.h>
#include "memory_ns.h"
#include "output.h"
#include <string>

#ifdef H5_LIB
#include "h5_c3po.h"
using namespace H5_C3PO_NS;
#endif

using namespace C3PO_NS;
using namespace C3PO_MEMORY_NS;

runningStat::runningStat()
:
    size(1), //default size 
    count_(NULL),
    run_mean_(NULL),
    run_var_(NULL),
    num_(0),
    num_g(0)
{
   
   MPI_Comm_rank(MPI_COMM_WORLD,&me_);
   MPI_Comm_size(MPI_COMM_WORLD,&nprocs_);

   binCentersDumped_ = false;

} 

runningStat::~runningStat() 
{
    delete count_;
    delete run_mean_;
    delete run_var_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::allocateMem(int _size)
{
    size = _size;
    create< int  >(count_,    size);
    create<double>(run_mean_, size);
    create<double>(run_var_,  size);
    create<double>(variance_, size);

    clear();

//    printf("running stats: mem allocated! \n");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::initialize(int size, const char* dumpFormat_,  const char * directory, const char* OpName)
{
     dumpFormat = dumpFormat_;
     OpName_ = OpName;
     file_.assign(directory);
     allocateMem(size);
   
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::computeGlobals(bool overwrite)
{
 int    bufCount;
 double bufMean;
 double bufVar;
 double* var_=variance();
 std::string f(file_);
 
 //MPI_Barrier(MPI_COMM_WORLD);

#ifdef H5_LIB
if(!dumpFormat.compare("hdf5"))
{
   
  double globalStat[size][3];

  
  if(me_==0) 
  {
     if(overwrite)    
     {
      f.append("_global.h5");
      createH5file(f);
     }
     else
     {
      char buf[40];
      sprintf(buf,"_global_%i.h5",num_g);
      f.append(buf);
      createH5file(f);
     }
  }
   //MPI_Barrier(MPI_COMM_WORLD);
   for(int i=0;i<size;i++)
    { 
      globalStat[i][0]=0.;
      globalStat[i][1]=0.;
      globalStat[i][2]=0.;
      
      for(int j=0;j<nprocs_;j++)
      {
       if(me_==j)
        MPI_Send(&count_[i],1, MPI_INT, 0,0, MPI_COMM_WORLD);
       if (me_==0)       
        MPI_Recv(&bufCount, 1, MPI_INT, j,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       if(me_==j)
        MPI_Send(&run_mean_[i],1, MPI_DOUBLE, 0,0, MPI_COMM_WORLD);
       if (me_==0)       
        MPI_Recv(&bufMean, 1, MPI_DOUBLE, j,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       if(me_==j)
        MPI_Send(&var_[i],1, MPI_DOUBLE, 0,0, MPI_COMM_WORLD);     
       if (me_==0)       
        MPI_Recv(&bufVar, 1, MPI_DOUBLE, j,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

       if(me_==0) 
       {
        globalStat[i][0]+=bufCount;                                           // Ctot= C1 + C2 + C3...
        globalStat[i][1]+=bufCount*bufMean;                                    // Ctot*Mtot= C1*M1 + C2*M2 + C3*M3 +....
        
        
        globalStat[i][2]+=bufCount*(bufVar + bufMean*bufMean);     // Ctot*(Mtot^2 + Vtot) = C1(V1 + M1^2) + C2(... 
    
       }
      }
     
    if(me_==0) 
    {
     if(globalStat[i][0]>0)
     {
      globalStat[i][1]=globalStat[i][1]/globalStat[i][0];
      globalStat[i][2]=globalStat[i][2]/globalStat[i][0] - globalStat[i][1]*globalStat[i][1];
     }
    } 
   }
  if(me_==0)
  ThreeArrayToH5(f, OpName_, globalStat, size ); 
    
}  //end HDF5
 
#endif

if(!dumpFormat.compare("json"))
{
  double globalStat[3][size];
   
  std::vector<double*>      datavec;
  std::vector<std::string>  namevec;
  
  if(me_==0) 
    f.append("_global.json");
   
 
   for(int i=0;i<size;i++)
    { 
      globalStat[0][i]=0.;
      globalStat[1][i]=0.;
      globalStat[2][i]=0.;
      
      for(int j=0;j<nprocs_;j++)
      {
       if(me_==j)
        MPI_Send(&count_[i],1, MPI_INT, 0,0, MPI_COMM_WORLD);
       if (me_==0)       
        MPI_Recv(&bufCount, 1, MPI_INT, j,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       if(me_==j)
        MPI_Send(&run_mean_[i],1, MPI_DOUBLE, 0,0, MPI_COMM_WORLD);
       if (me_==0)       
        MPI_Recv(&bufMean, 1, MPI_DOUBLE, j,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       if(me_==j)
        MPI_Send(&var_[i],1, MPI_DOUBLE, 0,0, MPI_COMM_WORLD);     
       if (me_==0)       
        MPI_Recv(&bufVar, 1, MPI_DOUBLE, j,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

       if(me_==0) 
       {
        globalStat[0][i]+=bufCount;                               // Ctot= C1 + C2 + C3...
        globalStat[1][i]+=bufCount*bufMean;                       // Ctot*Mtot= C1*M1 + C2*M2 + C3*M3 +....
        
        
        globalStat[2][i]+=bufCount*(bufVar + bufMean*bufMean);    // Ctot*(Mtot^2 + Vtot) = C1(V1 + M1^2) + C2(... 
    
       }
      }
     
    if(me_==0) 
    {
     if(globalStat[0][i]>0)
     {
      globalStat[1][i]=globalStat[1][i]/globalStat[0][i];
      globalStat[2][i]=globalStat[2][i]/globalStat[0][i] - globalStat[1][i]*globalStat[1][i];
     }
   }
  }
  if(me_==0) 
  {
   datavec.push_back(globalStat[0]);
   namevec.push_back("count");
   datavec.push_back(globalStat[1]);
   namevec.push_back("mean");
   datavec.push_back(globalStat[2]);
   namevec.push_back("variance");
   
   Output::createQJsonArrays(f,OpName_, namevec,datavec, size ,overwrite); 
  }   
} //end JSON

 num_g++;
 if(!overwrite)
  clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::dumpBinCenters(double binLow, double delta)
{

    if( me_>0 || binCentersDumped_ )
        return;

    std::string f(file_);

    #ifdef H5_LIB
    if(!dumpFormat.compare("hdf5"))
    {
        double centers[size];
        f.append("_binCenters.h5");
        createH5file(f);

        for(int i=0;i<size;i++)
        { 
            centers[i] = binLow + (double(i)+0.5) * delta;
        }

        OneArrayToH5(f, OpName_, centers, size ); 
    
    }  //end HDF5
    #endif

    if(!dumpFormat.compare("json"))
    {
        double centers[size];
        std::vector<double*>      datavec;
        std::vector<std::string>  namevec;
        f.append("_binCenters.json");
   
        for(int i=0;i<size;i++)
        { 
            centers[i] = binLow + (double(i)+0.5) * delta;
        }

        datavec.push_back(centers);
        namevec.push_back("centers");

        Output::createQJsonArrays(f,OpName_, namevec,datavec, true, size); 
    }   //end JSON

    binCentersDumped_ = true;

} 

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::add(double x, int bin)
{
    if(bin>=size)
     return;
    count_[bin]    += 1;
    double old      = run_mean_[bin];
    run_mean_[bin] += (x - run_mean_[bin])/count_[bin];
    run_var_[bin]  += (x - old)*(x-run_mean_[bin]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::updateFiles(bool overwrite)
{
 
  double* var_=variance();
  double data_[size][3];
   
  std::vector<double*>      datavec;
  std::vector<std::string>  namevec;
  
  std::string f(file_);
  
   #ifdef H5_LIB
    if(!dumpFormat.compare("hdf5"))
    {  
     if(overwrite)    
     {
      f.append(".h5");
      createH5file(f);
     }
     else
     {
      char buf[40];
      sprintf(buf,"_%i.h5",num_g);
      f.append(buf);
      createH5file(f);
     }
    }
   #endif
   
   if(!dumpFormat.compare("json"))
    f.append(".json");

   for(int i=0;i<size;i++)
    { 
      data_[i][0]=count_[i];
      data_[i][1]=run_mean_[i];
      data_[i][2]=var_[i];
      
      if(!dumpFormat.compare("json"))
      {
       datavec.push_back(data_[i]);
       namevec.push_back("bin");
      }
    }

   #ifdef H5_LIB
   if(!dumpFormat.compare("hdf5"))
    ThreeArrayToH5(f, OpName_, data_, size );
   #endif
    if(!dumpFormat.compare("json"))
     Output::createQJsonArrays(f,OpName_, namevec,datavec,overwrite,size); 
     
 num_++;
  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::clear()
{
   for(int i=0;i<size;i++)
   {
        count_[i]    = 0;
        run_mean_[i] = 0.0;
        run_var_[i]  = 0.0;
   }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
int* runningStat::count() 
{

    return count_;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * *    
double* runningStat::mean()
{
    return run_mean_;
}    

// * * * * * * * * * * * * * * * * * * * * * * * * * *
double* runningStat::variance()
{

    for(int i=0;i<size;i++)
        variance_[i] =  (count_[i] > 1) ? run_var_[i]/(count_[i]-1) : 0.0 ;
    
    return variance_; 
}


