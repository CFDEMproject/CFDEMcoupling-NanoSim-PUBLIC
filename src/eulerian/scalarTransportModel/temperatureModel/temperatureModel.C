/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "temperatureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(temperatureModel, 0);

addToRunTimeSelectionTable
(
	scalarTransportModel,
	temperatureModel,
	dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
temperatureModel::temperatureModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    scalarTransportModel(dict,sm),
    T_
    (
        IOobject
        (
            "T",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    Tsource_
    (
        IOobject
        (
            "Tsource",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    DT_(dict.lookup("DT"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

temperatureModel::~temperatureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void temperatureModel::createFields()
{}

void temperatureModel::update()
{
    //==============================
    // get references
    const volScalarField& voidfraction(particleCloud_.mesh().lookupObject<volScalarField> ("voidfraction"));
    const surfaceScalarField& phi(particleCloud_.mesh().lookupObject<surfaceScalarField> ("phi"));
    //==============================

    particleCloud_.forceM(1).manipulateScalarField(Tsource_,-1); // update source field 

    // solve scalar transport equation
    fvScalarMatrix TEqn
    (
       fvm::ddt(voidfraction,T_) - fvc::ddt(voidfraction,T_) //fvm::Sp(fvc::ddt(voidfraction,T_)) // PROBLEM: fvm::Sp does not work on volScalarField& ???? 
     + fvm::div(phi, T_) - fvm::Sp(fvc::div(phi),T_)
     - fvm::laplacian(DT_*voidfraction, T_)
     ==
       Tsource_
    );
    TEqn.relax();
    TEqn.solve();
}

volScalarField& temperatureModel::sourceField()
{
    return Tsource_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
