Sampling command specifiers
======================


_Sampling command specifiers_always need to be followed by the _command name_. Any _operation_ declared using a _sampling command specifier_ is called _sampling operation_ and has to be defined in the `c3po.json` file.

Syntax  
-------
```
...

 "sample2PC": {
                  "type"        : "TwoPointsCorr",
                  "marker1"     : "interFace",
                  "fieldToSample": "U_Favre",
                  "sampleCount" : 20,
                  "sampleDelta" : 0.5,
                  "save2Disk"   : true,
                  "save2Bin"    : true,
                  "lagrangian"  : false
              },
...
```
Every _sampling operation_ requires the user to specify some fields in the `c3po.json` file according to the "type" entry. However, any _sampling operation_ requires the following fields (for more details on type-specific entries see the corresponding *.md file):

* "type": requires a string of characters. The user can define which type of _sampling operation_ to perform.

* "save2Disk": requires a boolean value. If set to true, the sampled values will be written to disk in hdf5 or json format according to the specification provided in the "mainSettings" section of the `c3po.json` file. The corresponding output will consist of two arrays: the first is a marker, and the second one is the corresponding sampled value.

* "save2Bin": requires a boolean value. If set to true, the sampled values will be sent to the corresponding _binning operation_. If you choose this option, please be sure that a corresponding _binning operation_ exists. Otherwise, an error will be thrown.

Note, that documentation on how to use the different sampling operations delivered with CPPPO is provided in separate files.

Example
-------
```
...
                 
   "sample2PC": {
                  "type"        : "TwoPointsCorr",
                  "marker1"     : "interFace",
                  "fieldToSample": "U_Favre",
                  "sampleCount" : 20,
                  "sampleDelta" : 0.5,
                  "save2Disk"   : true,
                  "save2Bin"    : true,
                  "lagrangian"  : false
               },
  
   "samplePar": {
                  "type": "particles",
                  "save2Disk": true,
                  "save2Bin": true,
                 },
...
```

Read next
-----------
 - [general sampling](12_sampling_1_general.md)
 - [particle sampling](12_sampling_2_particle.md)
 - [two-points correlation sampling](12_sampling_3_twoPointCorr.md)
 - [angular sampling](12_sampling_4_angleVecVec.md)
 - [using the formula parser](12_sampling_5_formula.md)
 
Go back
-----------
 - [commands](10_commandTypes.md) 
