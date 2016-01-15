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

#include "CSVfieldOperations.h"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   fieldOperations Constructors
------------------------------------------------------------------------- */
CSVfieldOperations::CSVfieldOperations( CSVmesh  *        meshPtr, 
                                        double  **&     fieldsPtr
                                      )
:
fields_(fieldsPtr),
mesh_(meshPtr),
initialized_(true)
{
}

/* ----------------------------------------------------------------------
   fieldOperations Destructors
------------------------------------------------------------------------- */
CSVfieldOperations::~CSVfieldOperations()
{
}

/* ----------------------------------------------------------------------
   locating functions (IJK only)
------------------------------------------------------------------------- */
void CSVfieldOperations::cellIDtoIJK(int cellID_, int * IJK )
{
 
 int x=mesh_->NofCells()[0];    
 int y=mesh_->NofCells()[1];
 IJK[2] = int(cellID_ / (x*y)); //estimation of z-layer
 IJK[1] = int((cellID_-(x*y*IJK[2]))/x); //estimation of y-layer
 IJK[0] = int((cellID_-(x*y*IJK[2]))-((IJK[1]*x))); //estimation of x-layer

}

/*------------------------------------------------------------------------- */
int CSVfieldOperations::IJKtoCellID(int * IJK)
{
 
  int CellID_OF;          
   
    CellID_OF =   (IJK[2])*( mesh_->NofCells()[0] * mesh_->NofCells()[1] )
                + (IJK[1])* mesh_->NofCells()[0]
                +  IJK[0];
    
    return CellID_OF;

}

/*------------------------------------------------------------------------- */
void CSVfieldOperations::evaluateGradient( int      fieldID_ ,
                       double *     gradx,
                       double *     grady, 
                       double *     gradz
                     )
{
 //INPUT: 
 // fieldID_    ... id of the field to be sampled (field MUST be scalar)
 //OUTPUT: 
 // gradx       ... field containing the x-components of the gradient (length = length of field)
 
 double * baseField_ = fields_[fieldID_]; //This is an array containing the original field
 
 int    IJKCurrCell[3];
 int    IJKMinusNeighbor[3];
 int    IJKPlusNeighbor[3];
 int    indexMinusNeighbor;
 double valueMinusNeighbor;

 //Loop over all cell Ids

  for (int idC=0; idC<(mesh_->NofCells()[0]*mesh_->NofCells()[1]*mesh_->NofCells()[2]); idC++)
  {
    //compute the IJK index
    cellIDtoIJK(idC, &(IJKCurrCell[0]));

    //compute the index of the neighbor (in -x, and +x)
    IJKMinusNeighbor[0] = IJKCurrCell[0]-1;
    IJKMinusNeighbor[1] = IJKCurrCell[1];
    IJKMinusNeighbor[2] = IJKCurrCell[2];
    
    indexMinusNeighbor = IJKtoCellID(IJKMinusNeighbor);

    //pull out values using the index
    valueMinusNeighbor = baseField_[indexMinusNeighbor];

    //fill gradient fiels
    gradx[idC] = 999; //TODO


  }

}
