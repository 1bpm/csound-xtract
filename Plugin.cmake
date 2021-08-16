set(PLUGIN_NAME csxtract)
set(INCLUDES ${CSOUND_INCLUDE_DIRS})
set(LIBS "")

# Dependencies
find_package(Xtract)
check_deps(XTRACT_FOUND)
list(APPEND LIBS ${XTRACT_LIBRARIES})
list(APPEND INCLUDES ${XTRACT_INCLUDE_DIRS})

# Source files
set(CPPFILES src/opcodes.cpp src/dtw.cpp)
make_plugin(${PLUGIN_NAME} "${CPPFILES}" ${LIBS})
target_include_directories(${PLUGIN_NAME} PRIVATE ${INCLUDES})
