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

* "lagrangian":           requires a boolean value. If set to false, every cell in the domain will be filtered, otherwise only cells at particle centres will be filtered (Note that the "lagrangian": false algorithm is faster).

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
             "lagrangian": true
            },
              
...
```

Go back
-----------
 - [commands](10_commandTypes.md) 

