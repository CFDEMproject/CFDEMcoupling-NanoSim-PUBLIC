Sampling command specifiers
======================


"type": "particles"
----------------------
"Particles" _sampling operations_ do not require any additional entry. A "particles" _sampling operation_ will just send particle data (position, velocity, force) to a binning operation, or will dump them to disk in a tidy way. The user should consider that samples from "particles" are saved using the following markers:

* 0,1,2: for particle positions.

* 3,4,5: for particle velocities.

* 6,7,8: for forces acting on particles;

Example
-------
```
...

   "samplePar": {
                  "type": "particles",
                  "save2Disk": true,
                  "save2Bin": true,
                 },
...
```

Go back
-----------
 - [sampling](12_sampling.md)
