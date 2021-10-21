# PB-CT-DF-signal project
Shared PB-CT-DF-signal project source files
PbiDF.cpp file contains code for extraction of dark-field signal from in-line phase-contrast images

Here are the steps:
1) Padded experimental projections to the nearest powers of 2 in both directions.
2) TIE-hom retrieval from the padded projections.
3) Forward propagation of the object-plane complex amplitude based on the TIE-Hom reconstructed amplitude and phase.
4) Subtraction of the forward-propagated intensity obtained at step 3 from the original projection.
5) Addition of constant 1 to the result of the previous step (to remove any negative values).
6) Born-hom reconstruction from the result of the previous step.
7) Born-hom-iter reconstruction using the image from step 5 as I1 and the retrieved intensity from step 6 as I0.
The code reads inputs parameters from the text file "PbiDF.txt". Here is a list of inputs:
-	A path to the experimental projections
-	A path to output files (i.e. dark-field images)
-	Number of CT projections and stride
-	A range of the CT angle scan
-	Trim parameters to obtain projections with pixel numbers equal to the nearest powers of 2 in both directions.
-	Wavelength of the X-ray beam
-	Propagation distance (i.e. the distance from the sample to the detector)
-	Delta/beta values for both TIE-hom and Born approximations. These parameters are explained in the paper.  
