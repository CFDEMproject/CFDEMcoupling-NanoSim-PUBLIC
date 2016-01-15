Filtering command specifiers
======================


_Filtering command specifiers_always need to be followed by the _command name_. Any _operation_ declared using a _filtering command specifier_ is called _filtering operation_ and has to be defined in the c3po.json file.

Syntax  
-------
```
...

"favreU": 
            {
             "type":"Favre",
             "VectorfieldsToFilter": "U",
             "ScalarfieldsToFilter": " ",
             "phaseFractionField": "interFace",
             "lagrangian": false
            },
...

```
Every _filtering operation_ requires the user to fill some fields in the c3po.json file according to the "type" entry. However, any _filtering operation_ requires the following fields:

* "type":                 requires a string of characters. The user can define which type of _filtering operation_ to perform.

* "VectorfieldsToFilter": requires a string of characters. It has to be filled with the names of the vector fields to filter. 

* "ScalarfieldsToFilter": requires a string of characters. It has to be filled with the names of the scalar fields to filter. 

* "lagrangian":           requires a boolean value. If set to _false_, every cell in the domain will be filtered, otherwise only cells at particle/probe centres will be filtered (Note that the "lagrangian": false algorithm is faster). If this is set to _true_ the name of the probes/particles group should be provided with an additional entry:
 * "probesName": requires a string of characters. The user can has to enter the probes/particles group that will be used by this operation. Note that just one group per operation can be defined.

"type":"Favre" 
-------
_Filtering operations_ of type "Favre" will compute a phase-weighted average using the data specified in the `phaseFractionField` field (see below for details). Thus, these operations require an additional field to fill:

* "phaseFractionField": requires a string of characters. The user has to specify, in this field, the name of the the scalar field representing the phase fraction that is used for weighting the field to be averaged. Note, that the weighting field that was used to compute the Favre-averaged value will be also dumped to disk, and that it will be called according to the name specified as `phaseFractionField`. Of course, this field is NOT Favre-averaged. Note, in case `invertPhaseFraction` is specified as `true`, the (filtered) 1-`phaseFractionField` will be dumped.

If the user is using a two fluid model and wants to perform the "Favre" _filtering operation_ using a weighting field of the other phase, the user can use the following switch:

* "invertPhaseFraction": requires a boolean value. If set to `true`, 1-"phaseFractionField" will be used instead of "phaseFractionField" to weigth the field to be averaged. If not specified, the default value is `false`, i.e., the field specified as `phaseFractionField` will be used to weight the specified fields.

Note, that the setting for `invertPhaseFraction` will be applied to ALL fields within that filtering operation. Hence, the user has to specify at least two filtering operations in case he would like to perform a Favre-filtering on the solid and fluid phase.

"type":"Algebraic" 
-------
_Filtering operations_ of type "Algebraic" will calculate an arithmetic average. These kind of operations do not require any additional field to fill.


"type":"FavreRunningVariance" 
-------
_Filtering operations_ of type "FavreRunningVariance" will compute a phase-weighted average using the data specified in the `phaseFractionField` field (see below for details). Use this in case you also want to compute the variance.

Computing the variance of scalar, vector and mixed fields 
-------
Any filtering operation allows to compute variance fields. In order to do that, several new entries need to be provided:
* "VectorfieldsForVarianceName1" : requires a string of characters. This is the list of vector fields for variance calculation.
* "ScalarfieldsForVarianceName1" : requires a string of characters. This is the list of scalar fields for variance calculation.
* "VectorfieldsForVarianceComputeOffDiagonal" : requires an array of bools. If any element set to 'false', only the diagonal components of the variance field for the corresponding field in "VectorfieldsForVarianceName1" will be computed.
* "ScalarfieldsForVectorScalarMixedVariance" : requires a string of characters. This is the list of scalar fields for vector-scalar correlation. Enter "off" to disable.
* "VectorfieldsForVarianceName2" : requires a string of characters. This is the list of vector fields for vector-vector correlation.
* "ScalarfieldsForVarianceName2" : requires a string of characters. This is the list of scalar fields for scalar-scalar correlation.

Example
-------
```
...

"Favre": 
            {
             "type":"Favre",
             "VectorfieldsToFilter" : "U",
             "ScalarfieldsToFilter" : " ",
             "phaseFractionField"   : "voidage",
             "lagrangian"           : false,
             "invertPhaseFraction"  : false
            },
   
"Alg": 
            {
             "type":"Algebraic",
             "VectorfieldsToFilter": " ",
             "ScalarfieldsToFilter": "interFace",
             "lagrangian": true,
             "probesName":"samples"
            },
"Favre1": 
            {
             "type":"Favre",
             "VectorfieldsToFilter": "U",
             "ScalarfieldsToFilter": " ",
             "phaseFractionField": "voidfraction",
             "VectorfieldsForVarianceName1": "U",
             "VectorfieldsForVarianceComputeOffDiagonal": [false],
             "ScalarfieldsToFilter": " ",
             "ScalarfieldsForVectorScalarMixedVariance"  : " off",
             "ScalarfieldsForVarianceName1" : "",
             "ScalarfieldsForVarianceName2" : "",
             "lagrangian": true,
             "probesName":"particleCenter"
            }

...
```


Go back
-----------
 - [commands](10_commandTypes.md) 

Read more
----------
 - [probes and particles](20_probesAndParticles.md) 

