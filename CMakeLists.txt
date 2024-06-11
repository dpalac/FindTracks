# CMakeLists.txt for the "event" package. It creates a library and a main program.
# If ROOT is not installed in a default system location you need to tell CMake where to find it.
# Sourcing `thisroot.sh` already sets the required environment variables.
# Otherwise, you must tell the build system where to look for ROOT,
# for example by passing `-DROOT_DIR="/path/to/root/installation` at CMake configuration time.

   cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
   project(event)

# Locate the ROOT package and define a number of useful targets and variables.
   find_package(ROOT REQUIRED COMPONENTS RIO Net)
   find_package (Eigen3 3.4 REQUIRED NO_MODULE)
   if (TARGET Eigen3::Eigen)
    # Use the imported target
   endif (TARGET Eigen3::Eigen)
   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
   file(GLOB SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
   

   message(STATUS "> ${CMAKE_CURRENT_SOURCE_DIR}/inc")
   message(STATUS "> ${CMAKE_SOURCE_DIR}")

# Create a shared library.
# Passing cmake targets such as `ROOT::RIO` as dependencies (rather than plain
# library names for example via ${ROOT_LIBRARIES}) ensures that properties such as required
# include directories and C++ standard are propagated to our libraries or executables.
# Note: To ensure compatibility with Cling, targets *must* be compiled using the
# same C++ standard as ROOT was compiled with.
   add_library(RNSC SHARED "${SRCS}")

   target_link_libraries(RNSC PUBLIC Eigen3::Eigen ROOT::RIO ROOT::Net ROOT::Hist ROOT::Core ROOT::Imt ROOT::Graf ROOT::Graf3d ROOT::Gpad ROOT::ROOTVecOps ROOT::Tree ROOT::TreePlayer ROOT::Rint ROOT::Postscript ROOT::Matrix ROOT::Physics ROOT::MathCore ROOT::Thread ROOT::MultiProc ROOT::ROOTDataFrame)


# Create the main program using the library.
   add_executable(Tracks "${SRCS}")
   target_link_libraries(Tracks RNSC)
