# Fido

Code is in `fido/parallel_graph_align/main_code.cpp` and `fido/headers`.

Instructions to run the code: 

1. Download and install kokkos and kokkos-kernels
2. If using KokkosResilience fork then use tab 3.2.00 for kokkos-kernels
3. Create a build directory
4. Run CMake. ex
```
cmake \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DKokkos_DIR=/home/ntan1/KokkosResilience/kokkos/build/install/lib64/cmake/Kokkos \
  -DKokkosKernels_DIR=/home/ntan1/KokkosResilience/kokkos-kernels/build/install/lib64/cmake/KokkosKernels \
  ..
```
5. `make`
6. This will produce a binary called Fido
7. Run command: `mpirun -n N ./Fido path_to_graph1 path_to_graph2 ../fido/orbit_list.txt path_to_timing_file.txt --kokkos-threads=1` This will compare `graph1` and `graph2` using the orbits specified in `orbit_list.txt` and `N` processes.
8. Similarity matrix will be output to `out_similarity_matrix.txt`

