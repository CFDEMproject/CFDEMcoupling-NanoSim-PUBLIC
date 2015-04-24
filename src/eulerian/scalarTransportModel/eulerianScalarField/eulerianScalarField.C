/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2014- Stefan Radl, TU Graz, Austria

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "eulerianScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(eulerianScalarField, 0);

defineRunTimeSelectionTable(eulerianScalarField, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
eulerianScalarField::eulerianScalarField
(
    const dictionary&   dict,
    cfdemCloud&         sm,
    word                modelType,
    int                 modelID
)
:
    dict_(dict),
    particleCloud_(sm),
    fieldName_(modelType),
    m_
    (   IOobject
        (
            fieldName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    mSource_
    (   IOobject
        (
            fieldName_+"Source",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    ),
    fieldType_("undefined")
{


    if ( m_.dimensions() == dimensionSet(0, 0, 0, 1, 0) )
    {
        speciesID_ = -1;
        fieldType_ = "temperature";
        Info << "eulerianScalarField:: found a Temperature field! " << endl;
    }
    else
    {
        speciesID_ = modelID;
        fieldType_ = "species";
        Info << "eulerianScalarField:: found a species field, will assign speciesID: " 
             << speciesID_
             << endl;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

eulerianScalarField::~eulerianScalarField()
{}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void eulerianScalarField::update(surfaceScalarField phi, volScalarField voidfraction, volScalarField nuEff, scalar Sc) const 
{

    particleCloud_.forceM(1).manipulateScalarField(mSource_,speciesID_);

    // solve scalar transport equation
    fvScalarMatrix mEqn
    (
       fvm::ddt(voidfraction, m_) - fvc::ddt(voidfraction, m_)  
     + fvm::div(phi, m_)
     - fvm::Sp(fvc::div(phi), m_)
     - fvm::laplacian(nuEff/Sc*voidfraction, m_) 
     ==
       mSource_
    );
    mEqn.relax();
    mEqn.solve();


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
