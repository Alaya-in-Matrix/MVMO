cmake_minimum_required(VERSION 3.2.1)
project(mvmo)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
endif()
message(STATUS "MVMO Build type: ${CMAKE_BUILD_TYPE}")


# Suppress a warning from when compiling Eigen headers
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

set(MVMO_SRC MVMO.cpp MVMO_util.cpp)
set(MVMO_LIB mvmo)
add_library(${MVMO_LIB} ${MVMO_SRC})

set(DEMO_SRC demo.cpp)
set(DEMO_EXE demo)
add_executable(${DEMO_EXE} ${DEMO_SRC})
target_link_libraries(${DEMO_EXE} ${MVMO_LIB})

set_property(TARGET ${MVMO_LIB} ${DEMO_EXE} PROPERTY CXX_STANDARD 11)

# debug macro
option(MVMO_DEBUG "Debug macro for optimizer" OFF)
if(MVMO_DEBUG)
    message(STATUS "Open debug marcro")
    add_definitions(-DMYDEBUG)
    if(DEFINED RAND_SEED)
        message(STATUS "Rand seed is ${RAND_SEED}")
        add_definitions(-DDEBUG_RAND_SEED=${RAND_SEED})
    else()
        add_definitions(-DDEBUG_RAND_SEED=3)
    endif()
else()
    message(STATUS "Close debug marcro")
    remove_definitions(-DMYDEBUG)
    remove_definitions(-DDEBUG_RAND_SEED)
endif()

# if(NLOPT_PATH)
#     message(STATUS "NLOPT_PATH: ${NLOPT_PATH}")
#     include_directories("${NLOPT_PATH}/include")
#     find_library(NLOPT NAMES nlopt PATHS ${NLOPT_PATH}/lib NO_CMAKE_SYSTEM_PATH)
# else()
#     find_library(NLOPT NAMES nlopt)
# endif()
# if(NLOPT)
#     message(STATUS "Nlopt library: ${NLOPT}")
#     target_link_libraries(${MVMO_LIB} ${NLOPT})
# else()
#     message(FATAL_ERROR "NLOPT not found")
# endif()

# https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})


message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS ${MVMO_LIB} ${DEMO_EXE}
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
install(FILES MVMO.h MVMO_util.h
    DESTINATION inc 
    PERMISSIONS OWNER_READ GROUP_READ)
