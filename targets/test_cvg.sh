#!/bin/bash

## Perform the test on the circumcentric duals
rm -rf build
mkdir build
cd build
cmake ..
make -j12 test_square_cvg test_L_cvg test_tria_cvg
# Test the DEC convergence
echo "Testing Square ..."
./test_square_cvg ../meshes/square_centered.msh 8 > ../../dec_2d_square_centered_cvg_circ.txt
echo "Testing L-Domain ..."
./test_L_cvg ../meshes/L_centered.msh 8 > ../../dec_2d_L_centered_cvg_circ.txt
echo "Testing Triangle ..."
./test_tria_cvg ../meshes/tria.msh 8 > ../../dec_2d_tria_cvg.txt

cd ..
rm -rf build
