Specification of a 'scalarTransportModel'
=================

Syntax
-------------------------------

Defined in constant/scalarTransportProperties dictionary.

scalarTransportModel {theDesiredModel};


{theDesiredModel}Props
{
    {Provide here parameters for your model} 
};


Example 1: None - Do not use scalar transport
----------------------------

scalarTransportModel none;


Example 2: Multi-purpose transport model
----------------------------

scalarTransportModel generalManual;

generalManualProps
{
    ScT 0.7; //optional
    PrT 0.7; //optional

    //in case the user does not want to generate a separate field,
    //the volumetric heat capacity can be set here (as a global constant)
    cpVolumetric 1000; //this is the (mixture) density times the heat capacity
                       //must have dimensions [J/K/(m_voidspace)³]
                       //will only be used if cpVolumetricFieldName, or updateMixtureProperties = false
                       //see the input upon initialization of the transportModel

//    cpVolumetricFieldName   "cpRho";  //MUST NOT be set in case cpVolumetric is set
//    rhoMixFieldName         "rhoMix"; //Optional

    eulerianFields
    (
        C   //any concentration fields MUST be first to have correct numbering
        T
    );
    
    fvOptionsC {};  //Optional: any fvOptions for the transport model, must have correct name acc. to eulerianFields
    fvOptionsT {};  //Optional: any fvOptions for the transport model, must have correct name acc. to eulerianFields
};

Description
------------------------------------

In case a solver is used that relies on the'src/eulerian/scalarTransportModelsCFDEM/' library, the user MUST specify the keyword 'scalarTransportModel' in the constant/scalarTransportProperties files.
 

