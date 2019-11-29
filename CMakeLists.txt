cmake_minimum_required(VERSION 3.1) #Проверка версии CMake.

set(PROJECT_NAME 000_Lebedev_AllLab) #Создаем обычную локальную переменную с именем проекта
set(CMAKE_CXX_STANDARD 11)
project(${PROJECT_NAME}) # Название проекта

set(PROBLEM_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/_bin CACHE STRING "Output directory")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${PROBLEM_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROBLEM_OUTPUT_DIRECTORY})

set(GlobalSearch GlobalSearch) #Переменная с именем библиотеки
set(pugixml pugixml) #Переменная с именем библиотеки
option(USE_MPI "Set OFF to disable MPI" ON)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE) #Создается переменная CMAKE_CONFIGURATION_TYPES со значением "Debug;Release" и меткой/разделом Configs, переменная сохраниется в специальном файле доступном из всех файлов,
if(NOT CMAKE_BUILD_TYPE) # проверка имеющихся у компилятора вариантов сборки
  set(CMAKE_BUILD_TYPE Release) #Если нет то устанавливаем Release
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pugixml)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/)

if(WIN32)
    option(gtest_force_shared_crt "" TRUE)
endif(WIN32)

#include( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Configure.cmake" )

############################### GTest ###############################
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/gtest")
include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest/googletest/include")
#if( WIN32 )
    include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest-mpi/win")
#else( WIN32 )
    #include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest-mpi/linux")
#endif( WIN32 )

################################ MPI ################################

if( ${USE_MPI} )

    find_package( MPI )
    if( MPI_FOUND )
	message ("MPI FOUND!!!")
        include_directories( ${MPI_INCLUDE_PATH} )
    else( MPI_FOUND )
		message ("MPI Not Found")
        set( USE_MPI OFF )
    endif( MPI_FOUND )
endif( USE_MPI )

add_subdirectory(pugixml)
add_subdirectory(GlobalSearch)
add_subdirectory(MyFunction)
add_subdirectory(Rastrigin2)
add_subdirectory(GKLS)
add_subdirectory(Grishagin)
add_subdirectory(StronginC3)

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

add_executable(${PROJECT_NAME_STR} ${SRC_FILES})

#file(GLOB hdrs "*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
#file(GLOB srcs "*.cpp") #Добавляем в переменную srcs все файлы с расширением .cpp

#include_directories(lib1)
#add_executable(${PROJECT_NAME} ${srcs} ${hdrs}) #Создаем проект с именем из переменной ${PROJECT_NAME}, в него добавляются файлы из переменных ${srcs} и ${hdrs}


 #Собираем проект в поддиректории lib1
#target_link_libraries(${PROJECT_NAME} ${MyLib1}) #Подключаем к проекту с именем ${PROJECT_NAME} библиотеку с именем ${MyLib1}