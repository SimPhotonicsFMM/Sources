Copyright (c) 2023, Mondher Besbes (LCF / CNRS / IOGS)

SimPhotonics is a useful and powerful Matlab toolbox for the simulation of nanophotonic structures.
It is based on the Modal Fourier Method (well known as RCWA) which can be used 
to model multilayer structure, 1D and 2D periodic metamaterials.
SimPhotonics includes a specific feature that allows the design of nanoparticles 
by the use of the Gielis's SuperFormula.
This toolbox is the result of research works developed in the Charles Fabry laboratory 
and in particular the contributions of J.P. Hugonin.

Other Features:
--------------
+ Reflection and transmission of incident light (TE and TM calculation at the same time)
+ Dispersive materials (table array or handle function)
+ Reduced CPU time with y-axis symmetry (TE or TM polarization)
+ Quick calculation of a Bragg mirror
+ Intuitive field visualization
+ Parallel computation
+ Userfriendly examples

List of main functions: ("Sources" folder)
----------------------
- SetGeom :         define geometric parameters
- Spectrum :        compute the spectrum of reflectivity and transmittivity
- CalculFieldFMM :  calculate the field distribution
- VisuFieldFMM :    plot the distribution of the e.m. field 
- VisuMesh :        plot the discretization of the structure


List of some examples and tutorials ("Examples" folder)
-----------------------------------
* TutoSpectrum.mlx
* AntiReflectionCoating.mlx
* SPRBioSensorGrating1D and 2D

 
