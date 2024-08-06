cmake_minimum_required(VERSION 3.10)
enable_language(Fortran)
enable_language(CXX)
project(FTPYExample)

include(FortranCInterface)
include(CMakePrintHelpers)
FortranCInterface_VERIFY(CXX)


# Set C++ standard to at least C++11
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use an environment variable for the library installation path
# Fallback to /usr/local/lib if the environment variable is not set
set(LIB_INSTALL_DIR $ENV{CUSTOM_LIB_DIR})
if(NOT LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR "/usr/local/lib")
endif()

set(CMAKE_INSTALL_RPATH ${LIB_INSTALL_DIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Prefer the newer FindPython module for better handling of Python versions and libraries
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(yaml-cpp REQUIRED)

# Include directories for Python, Boost, and yaml-cpp
include_directories(${Python_INCLUDE_DIRS})
include_directories(${YAML_CPP_INCLUDE_DIR})



# Define the KamodoWrapper library target
add_library(KamodoWrapper SHARED KamodoWrapper.cpp KamodoWrapper.h)
add_library(FortranWrapper SHARED FortranWrapper.f90)

# Link libraries for KamodoWrapper, including yaml-cpp
target_link_libraries(KamodoWrapper ${Python_LIBRARIES} ${YAML_CPP_LIBRARIES}  -static-libstdc++ )
# Define the executable targets
add_executable(KamodoWrapperExample main.f90)

# Link executables with KamodoWrapper and the necessary libraries
target_link_libraries(KamodoWrapperExample FortranWrapper KamodoWrapper  ${Python_LIBRARIES} ${YAML_CPP_LIBRARIES}  ) 

cmake_print_variables(YAML_CPP_INCLUDE_DIR)
cmake_print_variables(YAML_CPP_LIBRARIES)
cmake_print_variables(Python_LIBRARIES)



# Installation rules
install(TARGETS KamodoWrapperExample
        RUNTIME DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/bin")


# Find the Python site-packages directory to install Python modules
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)




#
