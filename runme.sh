#!/bin/bash

## Perform the test on the circumcentric duals
rm -rf build
mkdir build
cd build
cmake ..
make -j4 test_square_cvg test_L_cvg test_tria_cvg
cd targets
# Test the DEC convergence
./test_square_cvg ../../meshes/square_centered.msh 8 > ../../dec_2d_square_centered_cvg.txt
./test_tria_cvg ../../meshes/tria.msh 10 > ../../dec_2d_tria_cvg.txt
#./test_L_cvg ../../meshes/L_centered.msh 8 > ../../dec_2d_L_centered_cvg_circ.txt


## Test the MG convergence
#./test_mg ../../meshes/square_centered.msh 8 > ../../decmg_2d_vcycle_square_centered_circ.txt
#./test_mg ../../meshes/L_centered.msh 8 > ../../decmg_2d_vcycle_L_centered_circ.txt
#./test_mg ../../meshes/tria.msh 8 > ../../decmg_2d_vcycle_tria.txt
#
### Perform the test on the barycentric duals
#cd ../..
#rm -rf build
#mkdir build
#cd build
#cmake -DBARYCENTRIC_DUAL=ON ..
#make -j12
#cd targets
#
#./test_square_cvg ../../meshes/square_centered.msh 8 > ../../dec_2d_square_centered_cvg_bary.txt
#./test_square_cvg ../../meshes/square_uncentered.msh 8 > ../../dec_2d_square_uncentered_cvg_bary.txt
#
#./test_L_cvg ../../meshes/L_centered.msh 8 > ../../dec_2d_L_centered_cvg_bary.txt
#./test_L_cvg ../../meshes/L_uncentered.msh 9 > ../../dec_2d_L_uncentered_cvg_bary.txt
#
#./test_mg ../../meshes/square_centered.msh 8 > ../../decmg_2d_vcycle_square_centered_bary.txt
#./test_mg ../../meshes/square_uncentered.msh 9 > ../../decmg_2d_vcycle_square_uncentered_bary.txt
#
#./test_mg ../../meshes/L_centered.msh 8 > ../../decmg_2d_vcycle_L_centered_bary.txt
#./test_mg ../../meshes/L_uncentered.msh 9 > ../../decmg_2d_vcycle_L_uncentered_bary.txt
#
#./test_mg ../../meshes/circle.msh 8 > ../../decmg_2d_vcycle_circle_bary.txt
#
#./test_mg ../../meshes/annulus.msh 5 5 > ../../decmg_2d_vcycle_annulus_bary.txt

## Clean up
cd ../..
rm -rf build

## Plot
python extract_errors.py dec_2d_square_centered_cvg 1
python extract_errors.py dec_2d_tria_cvg 2
