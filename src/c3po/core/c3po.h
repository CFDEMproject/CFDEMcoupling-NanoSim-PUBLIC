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
/*-----------------------------------------------------------------------------------
Description
	Core class for CPPPO. Provides all the functions needed to communicate with the
	interfaces and holds all tha main classes.
-----------------------------------------------------------------------------------*/


#ifndef C3PO_C3PO_H
#define C3PO_C3PO_H

#include "stdio.h"
#include "mpi.h"
#include <vector>

namespace C3PO_NS
{

class c3po
{
 friend class c3poBase;

 public:

  c3po(int narg, char **args, MPI_Comm communicator);
  ~c3po();

  void init();
  void report();
  
 private:

  void help();

  // class to register all members of this class

  class c3poBaseInterfaceVector *c3poBaseInterfaceVector_;

  // ***** CODE FRAMEWORK **************

  // Core Classes
  class Timer              *timer_;
  class Control       *control_;    // controller for stand-alone usage
  class Comm          *comm_;       // inter-processor communication
  class Input         *input_;      // input to simulation - input script processing
  class Output        *output_;     // simulation output
  class Error         *error_;      // error handling
  class OperationContainer *operationContainer_; // container for filtering/sampling/binning operations
  class SelectorContainer  *selectorContainer_;  // container for selection of cell and particle IDs
  class DataStorage   *dataStorageInternal_;   // storage for particle data, hosts several containers
  class c3poMesh            *mesh_;   // mesh data

  // ***** END CODE FRAMEWORK **************
    
  public:

    //Misc
    int selectorCount() const;
  
    //Top level run functions
    void   setCells(int * number_of_cells, int * global_number_of_cells, double * cellsize) const;

    void   runFilters(int id) const;

    void   runSampling() const;

    void   runBinning() const;
    
    void   pushA1DSample(int sampleId, double sampleValue, double x1) const;

    void   flush() const;  //general handle to flush/clear/clean all operation containers

    //TODO: return the instructions for sampling in a cleaner way
    int    sampleCount(int id) ;     
    double sampleDelta(int id) ; 

    //General settings
    bool   verbose() const;

  
    double meshCheckTolerance() const;
    double meshFilterWidthTolerance() const;
    bool   meshVerbose() const;
    
    void registerVF(std::string,double*,double*,double*,int sp=1);
    void registerSF(std::string,double*);
    void GlobalVF(std::string , double*,double*,double*, int sp=1);
    void GlobalSF(std::string, double* );
    void resetFields();
    void resetGlobalFields();
    
    std::string getVFnames(int);
    std::string getSFnames(int);
    
    int getVFnamesNumber();
    int getSFnamesNumber();
    
    int getVFnumber();
    int getSFnumber();
    
    double getMaxDomainfromJson(int i) const;
    double getMinDomainfromJson(int i) const;
    
    double cellSizefromJson(int i) const;
    
    
    void addParticle(double m,double* pos,double* vel,std::vector< double*>* force, double* torque = NULL);
    void deleteParticles();
    
    const char* getOpFilterName(int id); 
    void setTime(std::string t);
    
    int getOpFiltNum();
    int getOpSampNum();
    int getOpBinNum();
    
    void addFieldToConvert(std::string name_);
    void refreshFields();
    
    std::string getFilterName(int id);
    
    int numberOfFilters();
    
    std::vector<std::string> vectorFTF(int id);
    std::vector<std::string> scalarFTF(int id);

    std::vector<std::string> vectorFTFVariance(int id);
    std::vector<std::string> scalarFTFVariance(int id);
    
    void registerDomainInfo(double maxDomain[3],double minDomain[3],double maxDomainGlobal[3],double minDomainGlobal[3]);
    void registerCell(double cellV,const double* xcoord,const double* ycoord,const double* zcoord);
    
    void setNofCells(const int x) const;
    
    void preRunOperations() const;
    
    int getDomainDecomposition() const;
    
    bool writeFields() const;
    
    void csv_columnsToRead( int *columns) const;   
    
    void checkIJK(bool struct_)   const;
    
    void clearMesh() const;  
    
};

} //end c3po_NS

#endif

