add_library(MESH incidence.cpp bd.cpp mesh.cpp)
add_library(DEC_LIB dec.cpp cochain.cpp dirac.cpp)
add_library(MG hierarchy.cpp diracmg.cpp)
add_library(SPEC error.cpp)
add_library(BARY_LIB bary.cpp)

target_link_libraries(DEC_LIB MESH)
target_link_libraries(MG MESH DEC_LIB)
target_link_libraries(SPEC MESH DEC_LIB MG)
