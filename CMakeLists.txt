# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(Gem)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Build executable------------------------------------------------------------
add_executable(Gem gem.C)
target_link_libraries(Gem Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file MMTHGEM_dielectrics.dat MMTHGEM_fields.txt MMTHGEM_mesh.mphtxt)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
