cmake_minimum_required(VERSION 3.23.1)

# add these options before PROJECT keyword
# disallow in-source build:
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(gyrotests)

# compiler requirements
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -Wall")

file(GLOB app_files ${CMAKE_SOURCE_DIR}/apps/*.cc)
file(GLOB src_files ${CMAKE_SOURCE_DIR}/gyrotests/*.cc)

######################
# required libraries #
######################

# gyronimo requirement
set(gyrotests_includes ${gyrotests_includes} ${PROJECT_SOURCE_DIR})
set(gyronimo_dir /home/master/boris/gyrofork/gyronimo)
set(gyronimo_build /home/master/boris/gyrofork/build)
find_library(found gyronimo ${gyronimo_build})
  if(found)
    message(STATUS "Found ::gyronimo:: library")
  else()
    message(FATAL_ERROR "::gyronimo:: library not found, update gyronimo_dir path. See readme.md for installation details.")
  endif()
set(gyrotests_includes ${gyrotests_includes} ${gyronimo_dir} ${gyronimo_dir}/misc/apps/include)
set(gyrotests_depending_libs ${gyrotests_depending_libs} gyronimo)
link_directories(${gyronimo_build})

# gsl requirement
find_package(GSL REQUIRED)
  set(gyrotests_includes ${gyrotests_includes} ${GSL_INCLUDE_DIRS})
  set(gyrotests_depending_libs ${gyrotests_depending_libs} ${GSL_LIBRARIES})
  
  # boost requirement
find_package(Boost 1.73.0 REQUIRED)
  set(gyrotests_includes ${gyrotests_includes} ${Boost_INCLUDE_DIRS})
  set(gyrotests_depending_libs ${gyrotests_depending_libs} ${Boost_LIBRARIES})

# netcdf requirement
find_program(NCXX4_CONFIG "ncxx4-config")
  if(NCXX4_CONFIG)
    message(STATUS "Found ncxx4-config")
    execute_process(COMMAND "ncxx4-config" --libs OUTPUT_VARIABLE netcdf_cxx_libs)
    execute_process(COMMAND "ncxx4-config" --cflags OUTPUT_VARIABLE netcdf_cxx_flags)
  else()
    message(SEND_ERROR "ncxx4-config not found, install Unidata's netCDF-CXX4")
  endif()
  string(STRIP ${netcdf_cxx_libs} netcdf_cxx_libs)
  string(STRIP ${netcdf_cxx_flags} netcdf_cxx_flags)
  set(gyrotests_includes ${gyrotests_includes}  ${netcdf_cxx_flags})
  set(gyrotests_depending_libs ${gyrotests_depending_libs}  ${netcdf_cxx_libs})

add_library(gyrotests SHARED ${src_files})
target_include_directories(gyrotests PUBLIC ${gyrotests_includes})

message(${gyrotests_includes})
message(${gyrotests_depending_libs})

foreach(app_file ${app_files})
  get_filename_component(target ${app_file} NAME_WE)
  add_executable(${target} ${app_file})
  target_include_directories(${target} PUBLIC ${gyrotests_includes})
  target_link_libraries(${target} PRIVATE gyrotests ${gyrotests_depending_libs})
endforeach()

# add_executable(gyronimo-tests ${CMAKE_SOURCE_DIR}/app/gyronimo-tests.cc)
# target_include_directories(gyronimo-tests PUBLIC ${gyrotests_includes})
# target_link_libraries(gyronimo-tests PRIVATE gyrotests ${gyrotests_depending_libs})