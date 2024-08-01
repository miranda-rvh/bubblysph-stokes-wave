rm -rfv *.out

gfortran BUBBLE2VTU.f -o bub.out
gfortran GRID2VTU.f -o grid.out
gfortran SURF2VTU.f -o surf.out
gfortran SLICE2VTU.f -o slice.out
gfortran PART2VTU_2D.f -o part.out
