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

#include "generalManual.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(generalManual, 0);

addToRunTimeSelectionTable
(
	scalarTransportModel,
	generalManual,
	dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
generalManual::generalManual
(
    const dictionary& dict,
    cfdemCloud&       sm
)
:
    scalarTransportModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    eulerianFieldList_(propsDict_.lookup("eulerianFields")),
    ScT_(0.7),
    PrT_(0.7)
{

    propsDict_.readIfPresent("ScT", ScT_);
    propsDict_.readIfPresent("PrT", PrT_);
    
    Info << "Using ScT = " << ScT_ << " and PrT " << PrT_ << endl;

    eulerianFields_ = new autoPtr<eulerianScalarField>[eulerianFieldList_.size()];
    for (int i=0;i<eulerianFieldList_.size();i++)
    {
        eulerianFields_[i] = eulerianScalarField::New
        (
            propsDict_,
            sm,
            eulerianFieldList_[i],
            i
        );
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

generalManual::~generalManual()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void generalManual::createFields()
{}

void generalManual::update()
{
    //==============================
    // get references
    const surfaceScalarField& phi(particleCloud_.mesh().lookupObject<surfaceScalarField> ("phi"));
    const volScalarField& voidfraction(particleCloud_.mesh().lookupObject<volScalarField> ("voidfraction"));
    //==============================

    //Loop through all eulerian fields and update them
    for (int i=0;i<eulerianFieldList_.size();i++)
    {
            if(eulerianScalarF(i).fieldType()=="Temperature")
                eulerianScalarF(i).update(phi, voidfraction, particleCloud_.turbulence().nuEff(), PrT_);
            else 
                eulerianScalarF(i).update(phi, voidfraction, particleCloud_.turbulence().nuEff(), ScT_);
    }
}

// ************************************************************
volScalarField& generalManual::sourceField(int i)
{
    return eulerianScalarF(i).mSource();
}

// ************************************************************
const eulerianScalarField& Foam::generalManual::eulerianScalarF(int i)
{
     return eulerianFields_[i];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
