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
#include "mpi.h"
#include "comm.h"
#include "timer.h"
#include "memory_ns.h"
#include "output.h"
#include "stdio.h"
#include <cstring>
#include <stdlib.h>

using namespace C3PO_MEMORY_NS;

/* ---------------------------------------------------------------------- */

Timer::Timer(c3po *ptr) : c3poBase(ptr)
{
  create(array,TIME_N);
  init();
  stamp();
  synchronized_start(TIME_OVERALL);

  timeDir_ = "c3po_timings";
  timeDirIsCreated_ = false;

}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
   destroy(array);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
  for (int i = 0; i < TIME_N; i++) array[i] = 0.0;

}

/* ---------------------------------------------------------------------- */

void Timer::stamp() const 
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  previous_time = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp(int which) const 
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  double current_time = MPI_Wtime();
  array[which] += current_time - previous_time;
  previous_time = current_time;
}

/* ---------------------------------------------------------------------- */

void Timer::synchronized_start(int which)
{
 // MPI_Barrier(MPI_COMM_WORLD);
  array[which] -= MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::synchronized_stop(int which)
{
  //MPI_Barrier(MPI_COMM_WORLD);
  double current_time = MPI_Wtime();
  array[which] += current_time;
}

/* ---------------------------------------------------------------------- */
double Timer::elapsed(int which)
{
  double current_time = MPI_Wtime();
  return (current_time - array[which]);
}

/* ---------------------------------------------------------------------- */
void Timer::dumpStats()
{

  char buf[50];
  sprintf(buf,"timing_proc%i.json", comm().me());
  std::string fileName(buf);
  std::vector<std::string> names;
  names.push_back("TIME_OVERALL");
  names.push_back("TIME_SELECTOR");
  names.push_back("TIME_SELECTOR_BC");
  names.push_back("TIME_FILTER_MPI");
  names.push_back("TIME_FILTER");
  names.push_back("TIME_FILTER_TOTAL");
  names.push_back("TIME_OUTPUT");
  names.push_back("TIME_SYNC");
  names.push_back("TIME_SAMPLING");
  names.push_back("TIME_BINNING");
  
  std::vector<double*> data;
  for (int i=0; i<TIME_N;i++)
      data.push_back(&(array[i]));
  
  output().generateDir(timeDir_, timeDirIsCreated_);
  std::string fullFileName(timeDir_+"/"+fileName);
  output().createQJsonArrays(fullFileName,"c3po_timings",names,data,1,true);

}
