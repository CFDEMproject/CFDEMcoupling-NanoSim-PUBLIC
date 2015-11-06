#!/bin/bash

#This script will run the stokesFilter application
./Allclean

blockMesh
setFields
#run in parallel only when all selectors are cellUnstruct. This because the decomposition produces an unstructured mesh.
#decomposePar -force
#mpirun -np 3 stokesFilter -parallel 
stokesFilter
cd octave
#substiture 'processor0' with 'processor1' in postproc.m when run in parallel
octave postproc.m
evince validation.eps
evince solution.eps

