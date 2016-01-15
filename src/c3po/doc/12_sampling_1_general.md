Sampling command specifiers
======================


"type": "general"
----------------------

"General" _sampling operations_ will draw samples over all the Eulerian domain. The user can specify both the sampled field (which can be a vector or a scalar field) and the the _marker_ field (which has to be a scalar field). When choosing a vector field to sample, the user has also to specify the field component. Additional entries for "general" _sampling operations_ are:

* "VFieldToSample": requires a string of characters. The user has to specify the vector fields he/she wants to sample. At this point, the user can also use previously filtered fields. To access those fields, it is necessary to use the following syntax: (i) the name of the original field, and (ii) "\_", and finally (iii)  _command name_ corresponding to that _filtering operation_. 

* "component": requires an integer value. This field is only required if one or more vector fields are sampled. If set to 0, the x component will be sampled. If set to 1, the y component will be sampled. If set to 2, the z component will be sampled.

* "SFieldToSample": requires a string of characters. The user has to specify the scalar fields he/she wants to sample. At this point, the user can also use previously filtered fields. To access those fields, it is necessary to use the following syntax: (i) the name of the original field, and (ii) "\_", and finally (iii)  _command name_ corresponding to that _filtering operation_. 

* "marker": requires a string of characters. The user can define which field data is used to mark each sample (to be used for calculating conditional averages). Markers can also be previously filtered data fields, which can be selected similarly to the "fieldToSample" specification.

* "sampleCount": requires an integer. The user can specify the number of cells to sample. If set to -1, the whole domain will be sampled.

Note that this _sampling operation_ requires as many _binning operations_ as the number of fields to be sampled ( +1 if a _formula_ is used) when "samve2Bin" is set to 'true'.
This sampling operation allows multimarking in order to calculate conditional expected values of the sampled fields in an arbitrary parameter space. Results are dumped to disk in "c3po_binning" as explained in "13_binning.md". In addition a multimarking _general sampling operation_ must be connected to the corresponding multimarking _binning operation_.  

It is also possible to sample inside a specific region of the domain adding the following entries:
* "selective": requires a bool. If set to _true_ just a section of the domain will be sampled. 
* "max": requires an array of double values. The maximum box size.
* "min": requires an array of double values. The minimum box size.
 

Example
-------
```
...
                 
   "sampleGeneral": {
                  "type"        : "general",
                  "marker"     : "interFace interFace_Favre",
                  "VFieldToSample": "U_Favre",
                  "SFieldToSample": " p ",
                  "component"   : 0,
                  "sampleCount" : -1,
                  "save2Disk"   : true,
                  "save2Bin"    : true
              },
...
```
Go back
-----------
 - [sampling](12_sampling.md)
