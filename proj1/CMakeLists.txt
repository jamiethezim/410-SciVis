cmake_minimum_required(VERSION 3.10)

PROJECT(proj1)
SET(VTK_DIR /Users/jamiezimmerman/projects/buildVTK)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj1 proj1.cxx)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj1 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj1 vtkHybrid)
endif()


