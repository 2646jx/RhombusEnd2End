cmake_minimum_required(VERSION 3.11)
project(Cheetah)
include(CMakeFindDependencyMacro)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_compile_options(-fPIC)
add_compile_options(-D_GLIBCXX_USE_CXX11_ABI=0)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

## Dependencies
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(Eigen3 3.3.3 REQUIRED)
find_package(SEAL REQUIRED)
find_package(emp-tool REQUIRED)
find_package(emp-ot REQUIRED)
find_package(OpenMP)
include_directories(${EIGEN3_INCLUDE_DIR}
                    ${EMP-OT_INCLUDE_DIRS} 
                    ${EMP-TOOL_INCLUDE_DIRS})

# seal
set(SEAL_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/include/SEAL-3.7)
set(SEAL_LIBRARIES ${CMAKE_BINARY_DIR}/deps/SEAL/lib/libseal-3.7.a)  

# SCI
set(SCI_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/SCI/src/)

## cheetah-core
add_library(gemini SHARED)
include_directories(include)
add_subdirectory(include)
target_link_libraries(gemini SEAL::seal Eigen3::Eigen Threads::Threads OpenSSL::Crypto)

if (USE_APPROX_RESHARE)
  target_compile_definitions(gemini PUBLIC USE_APPROX_RESHARE=1)
else()
  target_compile_definitions(gemini PUBLIC USE_APPROX_RESHARE=0)
endif()

## SCI
include_directories(Rhombus/)
include_directories(SCI/src/)
add_subdirectory(SCI)


## Rhombus
set(RHOMBUS_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/Rhombus")
set(RHOMBUS_LIBRARIES "${CMAKE_BINARY_DIR}/lib/libRhombus.a")
add_subdirectory(Rhombus)


## demo
macro (add_network_cheetah _name)
  add_executable(${_name}-cheetah "networks/main_${_name}.cpp")
  target_include_directories(${_name}-cheetah PUBLIC ${RHOMBUS_INCLUDE_DIRS})
  target_link_libraries(${_name}-cheetah PUBLIC gemini SCI-Cheetah Rhombus)
  target_compile_definitions(${_name}-cheetah PUBLIC SCI_OT=1 USE_CHEETAH=1)
endmacro()

macro(add_network_rhombus _name)
  add_executable(${_name}-rhombus "networks/main_${_name}.cpp")
  target_include_directories(${_name}-rhombus PUBLIC ${RHOMBUS_INCLUDE_DIRS})
  target_link_libraries(${_name}-rhombus PUBLIC gemini SCI-Cheetah Rhombus)
  target_compile_definitions(${_name}-rhombus PUBLIC SCI_OT=1 USE_CHEETAH=1 USE_RHOMBUS=1)
endmacro()



add_network_cheetah(resnet50)
add_network_rhombus(resnet50)