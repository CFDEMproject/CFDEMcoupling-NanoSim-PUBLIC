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

    This is a general scalar transport model that is able to manage
    multiple scalar fields in a list. The user has to specify manually
    (in a dictionary) which fields need to be managed

    Contributor:
    Stefan Radl, TU Graz

Class
    generalManual

SourceFiles
    generalManual.C

\*---------------------------------------------------------------------------*/

#ifndef generalManual_H
#define generalManual_H

#include "scalarTransportModel.H"
#include "forceModel.H"
#include "eulerianScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class generalManual Declaration
\*---------------------------------------------------------------------------*/

class generalManual
:
    public scalarTransportModel
{

protected:

    dictionary                      propsDict_;

    const wordList                  eulerianFieldList_;

    autoPtr<eulerianScalarField>*   eulerianFields_;

    scalar                          ScT_; // - Turbulent Schmidt number (optional, default = 0.7)

    scalar                          PrT_; // - Turbulent Prandtl number (optional, default = 0.7)

public:

    //- Runtime type information
    TypeName("generalManual");


    // Constructors

    //- Construct from components
    generalManual
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

    ~generalManual();


    // Member Functions

    const                eulerianScalarField& eulerianScalarF(int);

    void                 createFields();

    void                 update();

    volScalarField&      sourceField(int i);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
