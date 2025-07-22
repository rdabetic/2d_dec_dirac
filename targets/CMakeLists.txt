add_executable(test_square_cvg test_square_cvg.cpp)
add_executable(test_tria_cvg test_tria_cvg.cpp)
add_executable(test_L_cvg test_L_cvg.cpp)
add_executable(test_mg test_mg.cpp)
add_executable(test test.cpp)
add_executable(test_interp test_interp.cpp)

target_link_libraries(test_square_cvg 
                      mfem Eigen3::Eigen MESH
                      DEC_LIB MG)
target_link_libraries(test_mg 
                      mfem Eigen3::Eigen MESH
                      DEC_LIB MG SPEC)
target_link_libraries(test_L_cvg
                      mfem Eigen3::Eigen MESH
                      DEC_LIB MG)
target_link_libraries(test
                      mfem Eigen3::Eigen MESH
                      DEC_LIB MG)
target_link_libraries(test_tria_cvg
                      mfem Eigen3::Eigen MESH
                      DEC_LIB MG BARY_LIB)
target_link_libraries(test_interp
                      mfem Eigen3::Eigen MESH
                      DEC_LIB MG BARY_LIB)
