cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME_STR src)

project (${PROJECT_NAME_STR})

set(PROBLEM_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/_bin CACHE STRING "Output directory")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${PROBLEM_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROBLEM_OUTPUT_DIRECTORY})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pugixml)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/GKLS/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Grishagin/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Hill/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/sample_src/Shekel/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)
					
add_library(${PROJECT_NAME_STR} SHARED ${SRC_FILES})

if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(CONFIG_OUTPUT_DIR ${PROBLEM_OUTPUT_DIRECTORY})
  set_property(TARGET ${PROJECT_NAME_STR} PROPERTY FOLDER "problems")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CONFIG_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
  target_compile_options(${PROJECT_NAME_STR} PUBLIC -std=c++11)
endif()

target_link_libraries(${PROJECT_NAME_STR} ${pugixml} ${EXAMIN_MPFR_LIB} ${EXAMIN_MPIR_LIB} ${DISLIN_LIB})


