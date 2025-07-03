# bubblysph-stokes-wave

BubblySPH code instructions.
This should give you a simulation with the same conditions as the Basilisk Stokes wave test case.

To generate input files, we use a small code called gen2D:
1) Navigate to source/gen/
2) Compile the gen2D with the "make" command <- this only needs doing the first time, or if you modify the source/gen/datclass.F90
3) Run it, with "./gen2D"
4) When prompted, enter the resolution (ratio of wavelength to particle spacing). 100 is pretty coarse, 300 is reasonable, 600 is pretty high resolution.

You will now have the relevant input files "IPART" and "IBOUND" in the main directory.

To compile bubblysph, in the main directory, use the command "make slipwalls=1"

To run it, use the command "./bubblysph"

The code will generate a set of output files in "data_directory". These need to be converted and then they can be loaded into paraview for visualisation. To convert the files:
1) navigate to the directory "vtk_conv"
2) Run "sh quickbuild.sh" <- this will compile some small codes to do the conversion. This compilation is only required the first time.
3) Convert the output files by running "./part.out"

You will now have a set of ".vtu" files in the directory "paraview_files"
You can load these into Paraview for visualisation.
