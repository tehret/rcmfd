project(RCMFD)
CMAKE_MINIMUM_REQUIRED(VERSION 1.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-aliasing")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-write-strings")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops -g")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops -g -fsanitize=address -static-libasan")

if(BICUBIC)
    message("-- Using bicubic interpolation.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBICUBIC")
endif()

message(${CMAKE_CXX_FLAGS})

include_directories(inverse_compositional)
add_library(inverseCompositional SHARED
    inverse_compositional/matrix.cpp
    inverse_compositional/transformation.cpp
    inverse_compositional/bicubic_interpolation.cpp
    inverse_compositional/zoom.cpp
    inverse_compositional/mask.cpp
    inverse_compositional/inverse_compositional_algorithm.cpp
    )

add_executable(main main.cpp matching.cpp utils/drawing.cpp utils/numerics1.cpp utils/utils.cpp utils/flimage.cpp utils/filter.cpp utils/sift.cpp utils/iio.c)
TARGET_LINK_LIBRARIES(main png tiff jpeg inverseCompositional)
