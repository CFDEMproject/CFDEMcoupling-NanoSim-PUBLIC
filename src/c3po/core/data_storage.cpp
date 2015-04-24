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


#include "data_storage.h"
#include "input.h"
#include "comm.h"
#include "output.h"
#include "mesh.h"
#include "error.h"
#include "operation_container.h"
#include "selector_container.h"
#include <sstream>
#include <stdlib.h>
#include "timer.h"

#ifdef H5_LIB
 #include "h5_c3po.h"
 using namespace H5_C3PO_NS;
#endif


using namespace C3PO_NS;

/* ----------------------------------------------------------------------
  DataStorage Constructor
------------------------------------------------------------------------- */

DataStorage::DataStorage(c3po *ptr) :
   c3poBase(ptr),
      nproc(comm().nprocs()),
   timeName_("0"),
   nbody_(0),
   nbody_all_(0),
   isAllocated_(false),
   MaxNofPar_(0)  
{
    haveFieldDir_ = false;
    haveParticleDir_ = false;
} 

/* ---------------------------------------------------------------------------*/ 
DataStorage::~DataStorage()
{
    deleteRMA();
    deleteFields();
    deleteParticles();   
   
}

/////////////////////////////////////////////////////////////////////////////
                          // MEMBER functions
/////////////////////////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   Read all necessary global properties so all MPI procs have them
------------------------------------------------------------------------- */

void DataStorage::read()
{
    //Allocate or Set Data arrays
    allocateMe();
}

/* ----------------------------------------------------------------------
   Write all necessary global properties so all MPI procs have them
------------------------------------------------------------------------- */

void DataStorage::write()
{

    //TODO
    //Write particle information to file
//    OperationProperties op(OPERATION_OUTPUT,false,false,false);

    //output().write_screen_one("DataStorage is now writing some output");
//    data().write(op);
}


/* ----------------------------------------------------------------------
   Scatter global properties so all MPI procs have them
------------------------------------------------------------------------- */

void DataStorage::scatter()
{
    if(comm().nprocs() > 1)
        error().throw_error_one(FLERR,"TODO: implement daat scattering if required");

   //TODO: in CustomValueTracker
   //make send buffer
   //pack buffer
   //unpack buffer
}

/* ----------------------------------------------------------------------
   Scatter global properties so all MPI procs have them
------------------------------------------------------------------------- */

void DataStorage::parallelize()
{
    if(comm().nprocs() > 1)
        error().throw_error_one(FLERR,"TODO: delete per-particle data if not in my subdomain");

}

/* ----------------------------------------------------------------------
   Intitialization phase - done before tun
------------------------------------------------------------------------- */

void DataStorage::init()
{

    return;
}

void DataStorage::allocateMe() const
{

    printf("DataStorage: will allocate memory for %d operations in operationContainer(). \n", operationContainer().operationCount());

    output().write_screen_one("\nDataStorage allocation completed.");

}


/* ----------------------------------------------------------------------
   Add particles / cells to memory
------------------------------------------------------------------------- */
void DataStorage::addParticle(double d, double* pos, double* vel, std::vector < double* >* force, double* torque)
{
 Particle* par_ = new Particle();
 particles_.push_back(par_);
 
 par_->setradius(d);
 par_->setpos(pos);
 par_->setvel(vel);
 par_->setforce(force);
 
 if(torque != NULL)
  par_->settorque(torque);
}

/* ---------------------------------------------------------------------------*/
void DataStorage::deleteParticles()
{
 for (unsigned int i=0;i<particles_.size();i++)
  delete particles_[i];
 particles_.clear();
}
/* ---------------------------------------------------------------------------*/
void DataStorage::gatherParticleData() const
{
 int tmp_[comm().nprocs()];
 int localNPar_ = particles_.size();
 
 MPI_Allgather(&localNPar_,1,MPI_INT,tmp_,1,MPI_INT,MPI_COMM_WORLD);
 
 //find the max
 for(int i=0;i<comm().nprocs();i++)
  if(MaxNofPar_<tmp_[i]) MaxNofPar_=tmp_[i];
 
 
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::writeFields(std::string OpName)
{

  if(!input().storageWriteFields())
      return;

  setFileName(OpName);
  output().generateDir(dirNameField_, haveFieldDir_);

  timer().stamp();
   
  int NX=mesh().NofCells();
  int fVF_size= fVF_.size();
  int fSF_size= fSF_.size();
  
 #ifdef H5_LIB 
 if( input().dumpFormat().compare("hdf5")==0 )
 {

  createH5file(filename_);
 
  for (int n=0;n<fVF_size;n++)    
  {
   double data[NX][3];
    
    for (int j = 0; j < NX; j++) 
	 for (int i = 0; i < 3; i++) 
	  data[j][i]=*(fVF_[n]->value(i,j));       
   std::string name(fVF_[n]->name());
   
   ThreeArrayToH5(filename_, name.c_str() , data, NX);
   
  }
  
  for (int n=0;n<fSF_size;n++)    
  {
   double data[NX];
    
    for (int j = 0; j < NX; j++) 
      data[j]=fSF_[n]->value()[j];       
   std::string name(fSF_[n]->name());
   
    OneArrayToH5(filename_, name.c_str() , data, NX);
   
  }

 }
 #endif
 if(input().dumpFormat().compare("json")==0)
 {

   std::vector<double**>    dataVPointer_;
   std::vector<std::string> dataNV_;
   
   std::vector<double*>     dataS_;
   std::vector<std::string> dataNS_;
   
   //save vector data
   for (int n=0;n<fVF_size;n++)     
   {
         dataVPointer_.push_back(fVF_[n]->values());
         std::string name(fVF_[n]->name());
         dataNV_.push_back(name);
   }
   std::string vname_(filename_);
   vname_.append("Vectors.json");
   if(dataNV_.size()>0)
     output().createQJsonArrays(vname_,dirNameField_,dataNV_, 
                              dataVPointer_, NX, 3, 
                              fVF_[0]->space()[0], true);

   //save scalar data
   for (int n=0;n<fSF_size;n++)    
   {
        std::string name(fSF_[n]->name());
        dataS_.push_back(fSF_[n]->value());
        dataNS_.push_back(name);   
   }
   std::string sname_(filename_);
   sname_.append("Scalars.json");
   if(dataNS_.size()>0)
     output().createQJsonArrays(sname_,dirNameField_,dataNS_,dataS_,NX,true);
    
 }

 timer().stamp(TIME_OUTPUT);
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::setCellWin(MPI::Win* x) const  { cellwin_ = x;}


/* ---------------------------------------------------------------------------*/ 
void DataStorage::writeParticles()
{

  if(!input().storageWriteParticles())
      return;

  setFileName("particles");
  output().generateDir(dirNameParticle_, haveParticleDir_);

  timer().stamp();

  int npar=particles_.size();
  if(npar==0) return;
  cout << "\n number of particles: " << npar;
  std::string file_(dirNameParticle_+"/time");
  file_.append(NameChanges("particles"));
  

  
  #ifdef H5_LIB 
  if(input().dumpFormat().compare("hdf5")==0)
  { 
   file_.append(".h5");
   createH5file(file_); 
  
  std::vector< double > data; 
  for (int par=0; par<npar; par++)
   data.push_back(*particles_[par]->getradius());
 
   OneArrayToH5(file_, "radius", &data[0], npar);
   data.clear();
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getpos()[0]);
   
   OneArrayToH5(file_, "pos_x", &data[0], npar);
   data.clear();
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getpos()[1]);
   
   OneArrayToH5(file_, "pos_y", &data[0], npar);
   data.clear();
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getpos()[2]);
   
   OneArrayToH5(file_, "pos_z", &data[0], npar);
   data.clear();
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getvel()[0]);
  
   OneArrayToH5(file_, "vel_x", &data[0], npar);
   data.clear();
  
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getvel()[1]);
   
   OneArrayToH5(file_, "vel_y", &data[0], npar);
   data.clear();
  
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getvel()[2]);
   
   OneArrayToH5(file_, "vel_z", &data[0], npar);
   data.clear();
   
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->gettorque()[0]);
  
   OneArrayToH5(file_, "torque_theta", &data[0], npar);
   data.clear();
  
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->gettorque()[1]);
   
   OneArrayToH5(file_, "torque_phi", &data[0], npar);
   data.clear();
  

  int NofForces_=particles_[0]->getNofForces();
  for(int force=0;force<NofForces_;force++) 
  {
   char buf[40];
   sprintf(buf,"force_%i",force);
   std::string dsetName_(buf);
   for (int par=0; par<npar; par++)
    data.push_back((particles_[par]->getforce(force))[0]);
 
   dsetName_.append("_x");
   OneArrayToH5(file_, dsetName_.c_str() , &data[0], npar);
   data.clear();
   for (int par=0; par<npar; par++)
    data.push_back((particles_[par]->getforce(force))[1]);
  
   dsetName_.assign(buf);
   dsetName_.append("_y");
   OneArrayToH5(file_,dsetName_.c_str() , &data[0], npar);
   data.clear();
   dsetName_.assign(buf);
   dsetName_.append("_z");
   for (int par=0; par<npar; par++)
    data.push_back((particles_[par]->getforce(force))[2]);
  
   OneArrayToH5(file_,dsetName_.c_str() , &data[0], npar);
   data.clear();
  }  
 }
 #endif
 
 if(input().dumpFormat().compare("json")==0)
 {
   file_.append(".json");
   std::vector<double*> dataV_;
   std::vector<double>  data;
   std::vector<std::string> dataN_;
   std::vector<int> datanum_;
   
   
   for (int par=0; par<npar; par++)
   {
    char buf[8];
    sprintf(buf,"%i",par);
    std::string particle(buf); 
    
    dataV_.push_back(particles_[par]->getradius());
    dataV_.push_back(particles_[par]->getpos());
    dataV_.push_back(particles_[par]->getvel());
    dataV_.push_back(particles_[par]->gettorque());
    
    dataN_.push_back(particle + "_radius");
    datanum_.push_back(1);
    dataN_.push_back(particle + "_position");
    datanum_.push_back(3);
    dataN_.push_back(particle + "_velocity");
    datanum_.push_back(3);
    dataN_.push_back(particle + "_torque");
    datanum_.push_back(2);
    
    int NofForces_=particles_[0]->getNofForces();
    for(int force=0;force<NofForces_;force++) 
    {
     char buf[40];
     sprintf(buf,"%s_force_%i",particle.c_str(),force);
     std::string dn_(buf);
     dataN_.push_back(dn_);
     dataV_.push_back(particles_[par]->getforce(force));
     datanum_.push_back(3);
    }
 
  output().createQJsonArrays(file_,"particles",dataN_,dataV_,-1,true, &datanum_);
  }

 }

 timer().stamp(TIME_OUTPUT);
}
//---------------------------------------------------------------//
void DataStorage::writeParticleFields(std::string OpName)
{
   if(!input().storageWriteParticles())
      return;

  setFileName(OpName);
  output().generateDir(dirNameParticle_, haveParticleDir_);

  timer().stamp();

  int npar=particles_.size();
  if(npar==0) return;
  int fVF_size= fVF_.size();
  int fSF_size= fSF_.size();
  
  std::string file_(dirNameParticle_+"/time");
  file_.append(NameChanges(OpName));
 
  #ifdef H5_LIB  
 if(input().dumpFormat().compare("hdf5")==0)
 { 
  file_.append(".h5");
  createH5file(file_); 
  
  std::vector< double > data; 
 
  for(int i=0;i<fVF_size;i++)
   for(int j=0;j<3;j++)
   {
    for (int par=0; par<npar; par++)
     data.push_back(particles_[par]->filteredVector(i)[j]);
    
    char buf[40];
    sprintf(buf,"%s_%i",fVF_[i]->name().c_str(),j);
    std::string dataset(buf);
    OneArrayToH5(file_, dataset.c_str(), &data[0], npar);
    data.clear();
   }

   for(int i=0;i<fSF_size;i++)
   {
    for (int par=0; par<npar; par++)
     data.push_back(*(particles_[par]->filteredScalar(i)));
    
    OneArrayToH5(file_, fSF_[i]->name().c_str(), &data[0], npar);
    data.clear();
   }
   
  
 } 
 #endif 
 

if(input().dumpFormat().compare("json")==0)
 {
  //  error().throw_error_one(FLERR,"output of particle data in JSON format not yet supported");
  file_.append(".json");
  
  std::vector<double*> dataV_;
  std::vector<std::string> dataN_;
 
  for(int i=0;i<fVF_size;i++)
   for(int j=0;j<3;j++)
   {
    double * data = new double[npar];

    for (int par=0; par<npar; par++)
     data[par]=(particles_[par]->filteredVector(i)[j]);
    
    char buf[40];
    sprintf(buf,"%s_%i",fVF_[i]->name().c_str(),j);
    std::string dataset(buf);
    dataN_.push_back(buf);
    dataV_.push_back(data);
   }

   for(int i=0;i<fSF_size;i++)
   {
    double * data = new double[npar];
    for (int par=0; par<npar; par++)
     data[par]=*(particles_[par]->filteredScalar(i));
    
    dataN_.push_back(fSF_[i]->name());
    dataV_.push_back(data);
    
   }
   output().createQJsonArrays(file_,"fields at particle centers",dataN_,dataV_,npar,true);
   
   for(unsigned int i=0;i<dataV_.size();i++)
    delete dataV_[i];
    
   dataV_.clear();
 }

 timer().stamp(TIME_OUTPUT);
}
/* ---------------------------------------------------------------------------*/ 
void DataStorage::addfVF(std::string name,double* x,double* y,double* z, int sp)
{
 //create the field
  filteredVectorField *v_ = new filteredVectorField(name, x,y,z,sp);

 //reset values 
  for(int i=0; i < mesh().NofCells();i++)
   for(int j=0;j<3;j++)
    *v_->value(j,i)=0;
  
  fVF_.push_back(v_);
  
  for(unsigned int par=0;par<particles_.size();par++)
   particles_[par]->addVector();

}


/* ---------------------------------------------------------------------------*/ 
void DataStorage::addRMAvF(std::string name,double* x,double* y,double* z, int sp)
{
 // MPI::Win* *win= new MPI::Win*[3];
 // int totcells_=selectorContainer().NofCells();
  
  filteredVectorField *v_ = new filteredVectorField(name, x,y,z,sp);
//  filteredVectorField *fv_ = new filteredVectorField(name, totcells_);

  RMAvF_.push_back(v_);
  
 // win[0]=comm().goRMA(x);
 // win[1]=comm().goRMA(y);
 // win[2]=comm().goRMA(z);
  
 // winV_.push_back(win);
 // fVF_.push_back(fv_);*/
  
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::addfSF(std::string name,double* x)
{
 //create the field
  filteredScalarField *s_ = new filteredScalarField(name, x);
 //reset values
  for(int i=0; i < mesh().NofCells();i++)  
   s_->value()[i]=0;
  
  fSF_.push_back(s_);
  
  for(unsigned int par=0;par<particles_.size();par++)
   particles_[par]->addScalar();

  

}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::addRMAsF(std::string name,double* x)
{
  //int totcells_=selectorContainer().NofCells();
  
  filteredScalarField *s_ = new filteredScalarField(name, x);
 //filteredScalarField *fs_ = new filteredScalarField(name, totcells_);

  RMAsF_.push_back(s_);
 // MPI::Win* win = comm().goRMA(x);
 // winS_.push_back(win);
 
 
// fSF_.push_back(fs_);
  
}

/* ---------------------------------------------------------------------------*/ 
filteredVectorField* DataStorage::fVF(std::string name)
{
  int fVF_size= fVF_.size();
  
  for (int i=0;i<fVF_size;i++)
  {
  
    if (fVF_[i]->name()==name)
     return fVF_[i];
 
  }
  
 char buf[100];
 sprintf(buf,"Can not find a Vector field named: %s inside C3PO, have you registered it?", name.c_str());
 error().throw_error_all("data_storage.cpp",-1,buf);
 return NULL;

}

/* ---------------------------------------------------------------------------*/ 
filteredScalarField* DataStorage::fSF(std::string name)
{
  int fSF_size= fSF_.size();
  
  for (int i=0;i<fSF_size;i++)
  {
  
    if (fSF_[i]->name()==name)
     return fSF_[i];
    
    
  }
  
  char buf[100];
  sprintf(buf,"Can not find a Scalar field named: %s inside C3PO, have you registered it?", name.c_str());
  error().throw_error_all("data_storage.cpp",-1,buf);
  return NULL;

}

/* ---------------------------------------------------------------------------*/ 
MPI::Win** DataStorage::getWinV(std::string name)
{
  int RMAvF_size= RMAvF_.size();
  
  for (int i=0;i<RMAvF_size;i++)
  {
  
    if (RMAvF_[i]->name()==name)
     return winV_[i];
  }
  error().throw_error_all("data_storage.cpp",0,"I can not find a window...");
  return NULL;
}

/* ---------------------------------------------------------------------------*/ 
MPI::Win* DataStorage::getWinS(std::string name)
{
  int RMAsF_size= RMAsF_.size();
  
  for (int i=0;i<RMAsF_size;i++)
  {
  
    if (RMAsF_[i]->name()==name)
     return winS_[i];
  }
  error().throw_error_all("data_storage.cpp",0,"I can not find a window...");
  return NULL;
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::deleteFields()
{
 

 int fVF_size= fVF_.size();
 for (int i=0;i<fVF_size;i++)
   delete fVF_[i];
   
 fVF_.clear();
 
 int fSF_size= fSF_.size();
 for (int i=0;i<fSF_size;i++)
   delete fSF_[i];
   
 fSF_.clear(); 

}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::deleteRMA()
{
  int winV_size = winV_.size();
  for (int i=0;i<winV_size;i++)
   for(int n=0;n<3;n++)
    {//(winV_[i])[n]->Free();
     //delete winV_[i][n];
     }
   //winV_.clear();

  int winS_size= winS_.size();
  for (int i=0;i<winS_size;i++)
   {
     //winS_[i]->Free();
    //delete winS_[i];
   }

  winS_.clear();   
 

  for(unsigned int i=0; i<RMAvF_.size();i++)
   delete RMAvF_[i];
  
  RMAvF_.clear();
  

  for(unsigned int i=0; i<RMAsF_.size();i++)
   delete RMAsF_[i];
  
  RMAsF_.clear();

}

/* ---------------------------------------------------------------------------*/ 
std::string DataStorage::NameChanges(std::string OpName)
{
   std :: ostringstream temp[3];
   
   for (int i=0;i<3;i++)
    temp[i] << selectorContainer().filterWidth()[i];
   
   std :: ostringstream temp2;
   temp2 << comm().me();
   
   std::string result("_"+timeName_+"_"+OpName.c_str()+"_processor"+temp2.str());
 
   return result;
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::setFileName(std::string OpName)
{
 dirNameField_.assign(   "c3po_dataStorage_fields");
 dirNameParticle_.assign("c3po_dataStorage_particles");

 if(input().dumpFormat().compare("hdf5")==0)
    filename_.assign(dirNameField_+"/results"+NameChanges(OpName)+".h5");
 
 if(input().dumpFormat().compare("json")==0)
    filename_.assign(dirNameField_+"/results"+NameChanges(OpName));

}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::convertFilterToRMA(std::string name_)
{
 //delete if already existent 
 for(unsigned int i=0;i<RMAvF_.size();i++)
  if(name_.compare(RMAvF_[i]->name())==0 )
   {
    delete RMAvF_[i];
    RMAvF_.erase(RMAvF_.begin()+i);   
   }
   
 for(unsigned int i=0;i<RMAsF_.size();i++)
  if(name_.compare(RMAsF_[i]->name())==0 )
   {
     delete RMAsF_[i];
     RMAsF_.erase(RMAsF_.begin()+i);  
   }
 
 //Create RMA field
 for(unsigned int i=0;i<fVF_.size();i++)
  if(name_.compare(fVF_[i]->name())==0)
  {
   addRMAvF(name_,fVF_[i]->value(0,0),fVF_[i]->value(1,0),fVF_[i]->value(2,0), *fVF_[i]->space());
   return;
  }
 for(unsigned int i=0;i<fSF_.size();i++)
  if(name_.compare(fSF_[i]->name())==0)
  {
   addRMAsF(name_,fSF_[i]->value());
   return;
  }
  
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::refreshRMAfields()
{
 for(unsigned int i=0; i<fieldsToConvert_.size(); i++)
  convertFilterToRMA(fieldsToConvert_[i]);
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::addFieldToConvert(std::string field_)
{
 //Check is the field is already added to list
 for(unsigned int i=0;i<fieldsToConvert_.size(); i++)
  if(field_.compare(fieldsToConvert_[i])==0) return;
 
 fieldsToConvert_.push_back(field_);
}
