cmake_minimum_required(VERSION 3.24)
project(laser_modulation)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math -mavx2 -g")  # -fopt-info-vec-optimized -fopt-info-vec-missed 

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find Boost 
find_package(Boost REQUIRED)

# Find FFTW3 (double-precision)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

add_executable(laser_modulation 
        src/main.cpp
        src/EvanescentField.cpp
        src/EvanescentField.h
        src/HelperFunctions.cpp
        src/HelperFunctions.h
        src/WavePacket.cpp
        src/WavePacket.h
        src/FourierTransformation.cpp
        src/FourierTransformation.h
        src/Grid2D.h
        src/Spline2D.h
        src/WavePacketModulation.cpp
        src/WavePacketModulation.h
        src/BilinearInterpolatorRegular.h
        src/BicubicInterpolator.h
        src/Timer.h
        src/Timer.cpp
)

message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")

# Include directories 
target_include_directories(laser_modulation PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
)

target_link_libraries(laser_modulation PRIVATE
        PkgConfig::FFTW
        OpenMP::OpenMP_CXX
)
