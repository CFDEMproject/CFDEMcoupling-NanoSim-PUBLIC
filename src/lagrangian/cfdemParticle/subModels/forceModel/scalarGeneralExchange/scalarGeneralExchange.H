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

    two way general scalar exchange between DEM and CFD
    convective heat and species transfer model following
    Li and Mason (2000), A computational investigation of transient heat
    transfer in  pneumatic transport of granular particles, Pow.Tech 112

    WARNING:
    This model REQUIRES the 'generalManual' speciesTransportModel

Class
    scalarGeneralExchange

SourceFiles
    scalarGeneralExchange.C

\*---------------------------------------------------------------------------*/

#ifndef scalarGeneralExchange_H
#define scalarGeneralExchange_H

#include "forceModel.H"
#include "averagingModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scalarGeneralExchange Declaration
\*---------------------------------------------------------------------------*/

class scalarGeneralExchange
:
    public forceModel
{
private:

    dictionary          propsDict_;

    const IOdictionary  scalarTransportProperties_;

    dictionary          generalPropsDict_;

    word                voidfractionFieldName_;

    word                velFieldName_;

    word                tempFieldName_;               
    word                partTempName_;
    word                partHeatFluxName_;
    word                partHeatTransCoeffName_;
    word                partHeatFluidName_;

    mutable double      **partDat_;                   // Lagrangian array
    mutable double      **partDatFlux_;               // Lagrangian array
    mutable double      **partDatTransCoeff_;         // Lagrangian array
    mutable double      **partDatFluid_;              // Lagrangian array
    
    bool                validPartFlux_;               //indicator if found, based on heat, but also used for species
    bool                validPartTransCoeff_;         //indicator if found, based on heat, but also used for species
    bool                validPartFluid_;              //indicator if found, based on heat, but also used for species

    scalar              lambda_;                      // fluid thermal conductivity [W/(m*K)]
    scalar              Cp_;                          // specific heat capacity [ J/(kg*K)]


    //Species Word Lists
    wordList            speciesFieldNames_;           //List with Eulerian fields to exchange (handed over)

    const wordList      partSpeciesNames_;
    const wordList      partSpeciesFluxNames_;
    const wordList      partSpeciesTransCoeffNames_;
    const wordList      partSpeciesFluidNames_;

    scalarList          DMolecular_;                 

    //Scalar properties
    scalar              maxSource_;                   // max (limited) value of src field

    mutable scalar      scaleDia_;

    void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("scalarGeneralExchange");

    // Constructors

        //- Construct from components
        scalarGeneralExchange
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~scalarGeneralExchange();


    // Member Functions
        void setForce() const;

        void manipulateScalarField(volScalarField&, int) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
