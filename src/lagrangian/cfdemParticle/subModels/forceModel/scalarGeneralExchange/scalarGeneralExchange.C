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

#include "error.H"

#include "scalarGeneralExchange.H"
#include "addToRunTimeSelectionTable.H"
#include "dataExchangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(scalarGeneralExchange, 0);

addToRunTimeSelectionTable
(
    forceModel,
    scalarGeneralExchange,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
scalarGeneralExchange::scalarGeneralExchange
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    scalarTransportProperties_                  //this is clumsy, but effective
    (
        IOobject
        (
            "scalarTransportProperties",
            sm.mesh().time().constant(),
            sm.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    generalPropsDict_(scalarTransportProperties_.subDict("generalManualProps")),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),             //common names/data
    velFieldName_(propsDict_.lookup("velFieldName")),
    tempFieldName_(propsDict_.lookup("tempFieldName")),                             //temperature names/data
    partTempName_(propsDict_.lookup("partTempName")),
    partHeatFluxName_(propsDict_.lookupOrDefault<word>(      "partHeatFluxName", "na")),
    partHeatTransCoeffName_(propsDict_.lookupOrDefault<word>("partHeatTransCoeffName", "na")),
    partHeatFluidName_(propsDict_.lookupOrDefault<word>(     "partHeatFluidName", "na")),
    partDat_(NULL),
    partDatFlux_(NULL),
    partDatTransCoeff_(NULL),
    partDatFluid_(NULL),
    validPartFlux_(false),
    validPartTransCoeff_(false),
    validPartFluid_(false),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    Cp_(readScalar(propsDict_.lookup("Cp"))),
    speciesFieldNames_( generalPropsDict_.lookup("eulerianFields")), 
    partSpeciesNames_(propsDict_.lookup("partSpeciesNames")),                       
    partSpeciesFluxNames_(propsDict_.lookup("partSpeciesFluxNames")),
    partSpeciesTransCoeffNames_(propsDict_.lookup("partSpeciesTransCoeffNames")),
    partSpeciesFluidNames_(propsDict_.lookup("partSpeciesFluidNames")),
    DMolecular_(propsDict_.lookup("DMolecular")),
    maxSource_(1e30),
    scaleDia_(1.)
{
    allocateMyArrays();

    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if(partHeatFluxName_!="na")         validPartFlux_=true;
    if(partHeatTransCoeffName_!="na")   validPartTransCoeff_=true;  
    if(partHeatFluidName_!="na")        validPartFluid_=true;

    if( validPartTransCoeff_ && !validPartFluid_ )
        FatalError <<"Transfer coefficient set, but and fluid name missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    

    if( !validPartTransCoeff_ && validPartFluid_ )
        FatalError <<"Fluid name set, but transfer coefficient missing. Check your entries in the couplingProperties! \n" 
                   << abort(FatalError);    
    
    if(!validPartFlux_ && !(validPartTransCoeff_ && validPartFluid_) )
        FatalError <<"You must set a valid heat flux name, or a valid transfer coefficient and fluid name \n" 
                   << abort(FatalError);

    // init force sub model
    setForceSubModels(propsDict_);

    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    forceSubM(0).readSwitches();

    particleCloud_.checkCG(true);

    if (propsDict_.found("scale"))
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));


    Info << "scalarGeneralExchange found the following speciesFieldNames: " << speciesFieldNames_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

scalarGeneralExchange::~scalarGeneralExchange()
{
    delete partDat_;
    delete partDatFlux_;

    if(validPartTransCoeff_)
       delete partDatTransCoeff_;

    if(validPartFluid_)
        delete partDatFluid_;

}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void scalarGeneralExchange::allocateMyArrays() const
{
    // Heat - get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partDat_,initVal,1);  // field/initVal/with/lenghtFromLigghts
    particleCloud_.dataExchangeM().allocateArray(partDatFlux_,initVal,1);

    if(validPartTransCoeff_)
        particleCloud_.dataExchangeM().allocateArray(partDatTransCoeff_,initVal,1);    
    if(validPartFluid_)
        particleCloud_.dataExchangeM().allocateArray(partDatFluid_,initVal,1);
}
// * * * * * * * * * * * * * * * public Member Functions  * * * * * * * * * * * * * //

void scalarGeneralExchange::setForce() const
{
    // do nothing
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void scalarGeneralExchange::manipulateScalarField(volScalarField& EuField, int speciesID) const
{
    //Set the names of the exchange fields
    word    fieldName;
    word    partDatName;
    word    partFluxName;
    word    partTransCoeffName;
    word    partFluidName;
    scalar  transportParameter;
    if(speciesID<0) //have temperature
    {
        fieldName          = tempFieldName_;
        partDatName        = partTempName_;
        partFluxName       = partHeatFluxName_;
        partTransCoeffName = partHeatTransCoeffName_;
        partFluidName      = partHeatFluidName_;
        transportParameter = lambda_;
    }
    else
    {
        fieldName          = speciesFieldNames_[speciesID];
        partDatName        = partSpeciesNames_[speciesID];
        partFluxName       = partSpeciesFluxNames_[speciesID]; 
        partTransCoeffName = partSpeciesTransCoeffNames_[speciesID]; 
        partFluidName      = partSpeciesFluidNames_[speciesID]; 
        transportParameter = DMolecular_[speciesID];
    }

    //==============================
    // get references
    const volScalarField& voidfraction_(particleCloud_.mesh().lookupObject<volScalarField> (voidfractionFieldName_));    // ref to voidfraction field
    const volVectorField& U_(particleCloud_.mesh().lookupObject<volVectorField> (velFieldName_));
    const volScalarField& fluidScalarField_(particleCloud_.mesh().lookupObject<volScalarField> (fieldName));            // ref to scalar field
    //==============================

    if (scaleDia_ > 1)
        Info << typeName << " using scale = " << scaleDia_ << endl;
    else if (particleCloud_.cg() > 1){
        scaleDia_=particleCloud_.cg();
        Info << typeName << " using scale from liggghts cg = " << scaleDia_ << endl;
    }

    // realloc the arrays
    allocateMyArrays();

    // reset Scalar field
    EuField.internalField() = 0.0;

    // get particle data
    particleCloud_.dataExchangeM().getData(partDatName,"scalar-atom", partDat_);

    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    // calc La based heat flux
    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar fluidValue(0);
    label  cellI=0;
    vector Us(0,0,0);
    scalar dscaled(0);
    scalar nuf(0);
    scalar magUr(0);
    scalar As(0);
    scalar Rep(0);
    scalar Pr(0);
    scalar Nup(0);
    scalar n = 3.5; // model parameter
    scalar sDth(scaleDia_*scaleDia_*scaleDia_);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> fluidScalarFieldInterpolator_(fluidScalarField_);

    scalar h1(0);
    scalar h2(0);
    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        //if(particleCloud_.regionM().inRegion()[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                if(forceSubM(0).interpolation())
                {
	                position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    fluidValue = fluidScalarFieldInterpolator_.interpolate(position,cellI);
                }else
                {
					voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                    fluidValue = fluidScalarField_[cellI];
                }

                // calc relative velocity
                Us      = particleCloud_.velocity(index);
                magUr   = mag(Ufluid-Us);
                dscaled = 2*particleCloud_.radius(index)/scaleDia_;
                As      = dscaled*dscaled*M_PI*sDth;
                nuf     = nufField[cellI];
                Rep     = dscaled*magUr/nuf;
                if(speciesID<0) //have temperature
                    Pr      = max(SMALL,Cp_*nuf*rhoField[cellI]/transportParameter); 
                else
                    Pr      = max(SMALL,nuf/transportParameter); //This is Sc for species

                if (Rep < 200)
                {
                    Nup = 2+0.6*pow(voidfraction,n)*sqrt(Rep)*pow(Pr,0.33); //This is Sh for species
                }
                else if (Rep < 1500)
                {
                    h1=pow(voidfraction,n);
                    h2=pow(Pr,0.33);
                    Nup = 2+0.5*h1*sqrt(Rep)*h2+0.02*h1*pow(Rep,0.8)*h2;
                }
                else
                {
                    Nup = 2+0.000045*pow(voidfraction,n)*pow(Rep,1.8);
                }
                scalar alpha = transportParameter*Nup/(dscaled);

                // calc convective heat flux [W]
                scalar tmpPartFlux     = alpha * As * (fluidValue - partDat_[index][0]);
                partDatFlux_[index][0] = tmpPartFlux;

                if(validPartTransCoeff_)
                    partDatTransCoeff_[index][0] = alpha;

                if(validPartFluid_)
                    partDatFluid_[index][0]      = fluidValue;


                if(forceSubM(0).verbose() && index >=0 && index <2)
                {
                    Pout << "index    = " <<index << endl;
                    Pout << "partFlux = " << tmpPartFlux << endl;
                    Pout << "magUr = " << magUr << endl;
                    Pout << "As = " << As << endl;
                    Pout << "r = " << particleCloud_.radius(index) << endl;
                    Pout << "dscaled = " << dscaled << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "Pr/Sc = " << Pr << endl;
                    Pout << "Nup/Shp = " << Nup << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "partDat_[index][0] = " << partDat_[index][0] << endl  ;
                    Pout << "fluidValue = " << fluidValue << endl  ;
                }
            }
        //}
    }

    particleCloud_.averagingM().setScalarSum
    (
        EuField,
        partDatFlux_,
        particleCloud_.particleWeights(),
        NULL
    );

    // scale with -1/(Vcell*rho*Cp) to give the source for the temperature field
    if(speciesID<0) //have temperature
        EuField.internalField() /= -rhoField.internalField()*Cp_*EuField.mesh().V();
    else
        EuField.internalField() /= -EuField.mesh().V();

    // limit source term
    scalar EuFieldInCell;
    forAll(EuField,cellI)
    {
        EuFieldInCell = EuField[cellI];

        if(mag(EuFieldInCell) > maxSource_ )
        {
             EuField[cellI] = sign(EuFieldInCell) * maxSource_;
        }
    }

    if(speciesID<0) //have temperature
        Info << "total convective particle-fluid heat flux [W] (Eulerian) = " 
             << gSum(EuField*rhoField*Cp_*EuField.mesh().V()) 
             << endl;
    else
        Info << "speciesID: " << speciesID 
             << ": total convective particle-fluid species flux [kmol/s] (Eulerian) = " 
             << gSum(EuField*1.0*EuField.mesh().V()) 
             << endl;

    // give DEM data
    if(validPartFlux_)
        particleCloud_.dataExchangeM().giveData(partFluxName ,
                                                "scalar-atom", 
                                                partDatFlux_);
    if(validPartTransCoeff_)
        particleCloud_.dataExchangeM().giveData(partTransCoeffName,
                                                "scalar-atom", 
                                                partDatTransCoeff_);
                
    if(validPartFluid_)
        particleCloud_.dataExchangeM().giveData(partFluidName,
                                                "scalar-atom", 
                                                partDatFluid_);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
