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
#include "selector_container.h"
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>
#include "input.h"
#include "output.h"
#include "operation_sampling_I_particles.h"

#ifdef H5_LIB
#include "h5_c3po.h"
using namespace H5_C3PO_NS;
#endif

using namespace C3PO_NS;


SamplingParticles::SamplingParticles(c3po *ptr,const char *name) 
: 
OperationSampling(ptr,name)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
SamplingParticles::~SamplingParticles()
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

void SamplingParticles::begin_of_step()
{
  (this->*run)();
}

//-----------------------------------------------------------------------------

void SamplingParticles::process_input(QJsonObject jsonObj)
{
 run=&SamplingParticles::particles;
}

//-----------------------------------------------------------------------------

void SamplingParticles::particles()
{
 int numpar_=dataStorage().numOfParticles();
 
 if(save2Bin_)
 {
  for(int par=0;par<numpar_;par++)
  {
   for(int j=0;j<3;j++)
   {
    insertSample(double(dataStorage().getParticle(par)->getpos()[j]),j);
    insertSample(dataStorage().getParticle(par)->getvel()[j],3+j);
    insertSample(dataStorage().getParticle(par)->getTotalForce(j),6+j);
   } 
   for(int j=0;j<2;j++)
   {
    insertSample(dataStorage().getParticle(par)->gettorque()[j],9+j);
   }
  }
  flushToBin();
 }
 
 
 if(save2Disk_)
  dataStorage().writeParticles();
 
 clearSamples();
} 

