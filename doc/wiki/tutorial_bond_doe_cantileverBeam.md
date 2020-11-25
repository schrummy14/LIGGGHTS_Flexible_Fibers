# LIGGGHTS Flexible Fibers

[Home](Home)

[How to Install](how_to_install)

[Commands](commands)

[Tutorials](tutorial_main_page)

## Bond Package DOE Cantilever Beam Tutorial

### Simulation Description

This tutorial shows how to set up a Design of Experiments (DOE) to determine
the Young's modulus and bond coefficent of a fiber utilizing the cantilever
beam test. The DOE will alter the bond Young's modulus and the bond damping
coefficient to create linear relationships between the frequency of oscilation
to the bond Young's modulus and the global decay of the oscilation to the bond
damping coefficient.

The user first creates the DOE parameters by running the python code
"*makeDOE*". This will create a file that will be read in by the LIGGGHTS
input script. The user can modify this file to add/modify the number of
parameters and the levels in which the DOE will run.

Once the simulations are done, the user may run the *extractProperties* file to
create linear models for the variables.
