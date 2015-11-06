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

Application / Class
    eulerianScalarField
    
Description
    Base class for moment transport model & FVM-based solver

\*---------------------------------------------------------------------------*/

#ifndef eulerianScalarField_H
#define eulerianScalarField_H

#include "cfdemCloud.H"
#include "forceModel.H"
#include "fvOptionList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class eulerianScalarField Declaration
\*---------------------------------------------------------------------------*/

class eulerianScalarField
{

protected:

    // Protected data
        const dictionary&       dict_;

        cfdemCloud&             particleCloud_;

        const word              fieldName_;

        mutable volScalarField  m_;         // transported quantity

        mutable volScalarField  mSource_;   // source of transported quantity

        mutable word            fieldType_;

        mutable fv::optionList  fvOptions_;

        int                     speciesID_;


public:

    //- Runtime type information
    TypeName("eulerianScalarField");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            eulerianScalarField,
            dictionary,
            (
                const dictionary&   dict,
                cfdemCloud&         sm,
                word                modelType,
                int                 modelID
            ),
            (dict,sm,modelType,modelID)
        );


    // Constructors

        //- Construct from components
        eulerianScalarField
        (
            const dictionary&   dict,
            cfdemCloud&         sm,
            word                modelType,
            int                 modelID
        );


    // Destructor

        virtual ~eulerianScalarField();


    // Selector

        static autoPtr<eulerianScalarField> New
        (
            const dictionary&   dict,
            cfdemCloud&         sm,
            word                modelType,
            int                 modelID
        );


    // Member Functions
    void update(surfaceScalarField, volScalarField, volScalarField, scalar Sc) const ;

    void bound(autoPtr<eulerianScalarField>*) const {};

    // Access
    word fieldType() const                 { return fieldType_; };

    inline volScalarField& m()       const {return m_;}; //returns the transported quantity

    inline volScalarField& mSource() const {return mSource_;}; //returns the source

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
