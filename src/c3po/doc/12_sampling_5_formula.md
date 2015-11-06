Sampling formula
======================


Description
----------------------

CPPPO allows the user to perform basic algebraic operations on sampled fields and to register the result with the operation marker. A _formula_ is an optional entry for every _sampling operation_.

* "Formula": requires a string of characters. The user can type the formula she/he wants CPPPO to execute.

The following rules applies to user defined formulas:

* The user can specify a numerator and a denominator with multiple operations for each. The expression for the numerator and denominator MUST be in brackets and separated with the '%' symbol.
* CPPPO will fill the data computed via the formular into a sample field with the extension 'formula0'
* elements specified via the keyword 'VFieldsToSample' can be accessed with 'vec1', 'vec2', ..., 'vec9'. It is not possible to use more than 9 vectors (or scalars) in such a way
* elements specified via the keyword 'SFieldsToSample' can be accessed with 'scalar1', 'scalar2', ..., 'scalar9'.
* Operations in the numerator and the denominator are executed from left to right. Thus, the order of the operations defines the output, and brackets in the numerator and denominator will be NOT respected!
* Thus, the formula "(vec1 - vec2) * scalar1" will give a different result as "scalar1 * (vec1 - vec2)"  
* The keyword "saveOnlyFormula" can be specified in case a formula is specified. If set to true, only the result of the formular will be saved (to bin or disk), but NOT the sampled vector and scalar fields.
* The user can also directly refer to sampled fields using numbers from 0 to 9. CPPPO will order the fields considering vector fields first using the ordering provided by the user input.
* The only operations allowed are:
 * '*' : multiplication
 * '/' : division
 * '+' : sum
 * '-' : difference

Note that the computed value will be considered as a new sample and, if "save2Bin" is set to 'true', it will require an additional _binning operation_. The user may want to set "saveOnlyFormula" : true in order to avoid specifying an uncessary large amount of bins.

Example
-------
```
...
                 
   "sampleGeneral": {
                  "type"        : "general",
                  "marker1"     : "interFace",
                  "VFieldsToSample": "U_Favre U_Alg ",
                  "SFieldsToSample": "rho",
                  "component"   : 0,
                  "Formula"     : "(vec1 - vec2) % ( scalar1 ) ", 
                  "saveOnlyFormula" : true,
                  "sampleCount" : -1,
                  "save2Disk"   : true,
                  "save2Bin"    : true
              },
...
```
In this example the formula entered corresponds to:
```
 (U_Favre-U_Alg) / rho 
```

Go back
-----------
 - [sampling](12_sampling.md)
