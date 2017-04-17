# Extended-Source-Diffraction

MATLAB code to simulate the diffraction of a wave through an arbitrary aperture. The mathematical justification of the code is included in a seperate pdf file.

To use, set the variables in diffraction_main.m for the size of the source and the properties of it. Then the aperture is set in aperture_funct.m, currently it is set up for a frensel zone plate but this can be altered for any aperture desired. 

The code can take a long time to run if a large aperture is used with a large source size, so some care is needed when setting up the problem to ensure the computer's memory isn't filled.

