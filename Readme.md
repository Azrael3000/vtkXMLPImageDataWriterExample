## Why:

VTK can be a bit difficult to work with, so here is a small example that writes out a \*.vti file in parallel. Each processor only writes part of the data and the data contains one cell scalar and one point scalar. So I hope it covers the functionality to get you started with more complex examples.

Of particular importance for the vtkXMLP\*Writer is the availability of a MultiProcessController which is also demonstrated in this little piece of code.

## Compile:

    mpicxx vtk.cpp -Wall -I/usr/include/vtk-6.2 -lvtkCommonCore-6.2 -lvtkImagingCore-6.2 -lvtkCommonExecutionModel-6.2 -lvtkIOCore-6.2 -lvtkIOXML-6.2 -lvtkIOParallelXML-6.2 -lvtkCommonDataModel-6.2 -lvtkParallelCore-6.2 -lvtkParallelMPI-6.2 -o vtk

Depending on your vtk version and installation folder this might need some modification.

## Execute:

    ./vtk

or

    mpirun -np 2 ./vtk

## Contribute:

If somebody finds an error or wants to make this more user friendly by providing a CMake script or comments to the code feel free to submit a pull request.
