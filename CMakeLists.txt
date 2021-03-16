project("csound_xtract")

cmake_minimum_required(VERSION 3.8)

set(APIVERSION "6.0")

# Release or Debug
set(CMAKE_BUILD_TYPE "Release")

# force make to print the command lines
set(CMAKE_VERBOSE_MAKEFILE on)

# path to Csound cmake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# set compilation flags
set(CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS}  "-std=c++11 -fpermissive -fPIC -w -DUSE_DOUBLE -DB64BIT")

# options
option(USE_LIB64 "Set to on to set installation dir for libs to lib64" OFF)
option(USE_DOUBLE "Use doubles for audio calculations" ON)
option(CPP11 "c++11" ON)

set(BUILDING_CSOUND_PLUGINS ON)

# ----------------------------------------------

include(FindCsound)
include(FindXtract)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

# -----------------------------------------------

function(addflag flag flagname)
  check_c_compiler_flag(${flag} ${flagname})
  if (${flagname})
    # message(STATUS "Setting C flag ${flag}")    
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}" PARENT_SCOPE)
  endif()
  check_cxx_compiler_flag(${flag} CXX_${flagname})
  if (CXX_${flagname})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
  endif()
endfunction(addflag)


# set optimization flags
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
  add_definitions(-fvisibility=hidden)
  if(NATIVE)
	add_definitions(-march=native)
  endif()
  
  include(CheckCCompilerFlag)
  include(CheckCXXCompilerFlag)

  addflag(-msse HAS_SSE)
  addflag(-msse2 HAS_SSE2)
  addflag(-mfgpath=sse HAS_FPMATH_SSE)
  
endif()

if(MINGW)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mstackrealign")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mstackrealign")
endif()

addflag(-ftree-vectorize HAS_TREE_VECTORIZE)
addflag(-ffast-math HAS_FAST_MATH)
addflag(-fomit-frame-pointer HAS_OMIT_FRAME_POINTER)


# -------------------------------------------------------------------

set(CS_FRAMEWORK_DEST "~/Library/Frameworks")


if(USE_LIB64)
  set(LIBRARY_INSTALL_DIR "lib64")
  add_definitions("-DLIB64")
else()
  set(LIBRARY_INSTALL_DIR "lib")
endif()

message(STATUS "LIBRARY INSTALL DIR: ${LIBRARY_INSTALL_DIR}")

# -------------------------------------------------------------------


if(USE_DOUBLE)
  message(STATUS ">>> using doubles")
  
  if(APPLE)
    set(CSOUNDLIB "CsoundLib64")
    set(PLUGIN_INSTALL_DIR "${CS_FRAMEWORK_DEST}/${CSOUNDLIB}.framework/Versions/${APIVERSION}/Resources/Opcodes64")
  else()
    set(CSOUNDLIB "csound64")
    set(PLUGIN_INSTALL_DIR "${LIBRARY_INSTALL_DIR}/csound/plugins64-${APIVERSION}")
  endif()
else()
  message(STATUS ">>> not using doubles")
  if(APPLE)
    set(CSOUNDLIB "CsoundLib")
    set(PLUGIN_INSTALL_DIR "${CS_FRAMEWORK_DEST}/${CSOUNDLIB}.framework/Versions/${APIVERSION}/Resources/Opcodes")
  else()
    set(CSOUNDLIB "csound")
    set(PLUGIN_INSTALL_DIR "${LIBRARY_INSTALL_DIR}/csound/plugins-${APIVERSION}")
  endif()
endif()


# -------------------------------------------------------------------
  
# Csound opcode build
find_package(Csound)


set(BUILD_PLUGINS_DIR ${CMAKE_CURRENT_BINARY_DIR})

if(NOT CSOUND_FOUND)
  message(FATAL_ERROR "Csound installation not found")
endif()


if(NOT XTRACT_FOUND)
  message(FATAL_ERROR "libxtract could not be found")
endif()

set(CPPFILES src/opcodes.cpp src/dtw.cpp)


include_directories(${CSOUND_INCLUDE_DIRS})
link_libraries(${XTRACT_LIBRARIES})
include_directories(${XTRACT_INCLUDE_DIRS})
include_directories(include)


add_library(csxtract SHARED ${CPPFILES})

set_target_properties(csxtract PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${BUILD_PLUGINS_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${BUILD_PLUGINS_DIR})
 
install(TARGETS csxtract LIBRARY DESTINATION "${PLUGIN_INSTALL_DIR}" )

