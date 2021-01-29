project(RCMFD)
CMAKE_MINIMUM_REQUIRED(VERSION 1.0)
 
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-aliasing")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wno-write-strings")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops -g")
##SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -funroll-loops -g -fsanitize=address -static-libasan")
#
#
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -fPIC")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -march=native -fPIC")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")


# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else ()
    set (CMAKE_C_STANDARD 99)
endif ()

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

add_library(rcmfd SHARED
    libmatching.cpp 
    matching.cpp
    utils/drawing.cpp
    utils/numerics1.cpp
    utils/utils.cpp 
    utils/flimage.cpp 
    utils/filter.cpp 
    utils/sift.cpp 
    )

add_library(wrcmfd SHARED
    wrapper_rcmfd.c
    )

target_link_libraries(wrcmfd rcmfd)
set_property(TARGET wrcmfd PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(TARGET wrcmfd POST_BUILD
   COMMAND gcc -I${CMAKE_CURRENT_SOURCE_DIR} -E ${CMAKE_CURRENT_SOURCE_DIR}/all.h >all.h
   COMMENT "Creating all.h"
)

add_executable(main main.cpp matching.cpp utils/drawing.cpp utils/numerics1.cpp utils/utils.cpp utils/flimage.cpp utils/filter.cpp utils/sift.cpp utils/iio.c)
TARGET_LINK_LIBRARIES(main png tiff jpeg inverseCompositional)
