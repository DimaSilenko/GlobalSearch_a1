file(GLOB hdrs "*.h*") #��������� � ���������� hdrs ��� ����� � ����������� .h
file(GLOB srcs "*.cpp")#��������� � ���������� srcs ��� ����� � ����������� .h


include(../cmake/examin_utils.cmake)
include(../cmake/examin_mpi.cmake)
include(../cmake/examin_lib.cmake)
include(../cmake/examin_doxygen.cmake)

examin_initialize()

setup_mpir(EXAMIN_MPIR_LIB_DIR EXAMIN_MPIR_INCLUDES EXAMIN_MPIR_LIB)
setup_mpfr(EXAMIN_MPFR_LIB_DIR EXAMIN_MPFR_INCLUDES EXAMIN_MPFR_LIB)
setup_mpreal(EXAMIN_MPREAL_INCLUDES)

include_directories(
  ${EXAMIN_MPIR_INCLUDES}
  ${EXAMIN_MPFR_INCLUDES}
  ${EXAMIN_MPREAL_INCLUDES}
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/SRC)

link_directories(
  ${EXAMIN_MPIR_LIB_DIR}
  ${EXAMIN_MPFR_LIB_DIR}
)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()



add_executable(${GlobalSearch} ${srcs} ${hdrs})

if( ${USE_MPI} )
	message("Set MPI START")
    if( MPI_COMPILE_FLAGS )
        set_target_properties( ${GlobalSearch} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}" )
    endif( MPI_COMPILE_FLAGS )

    if( MPI_LINK_FLAGS )
        set_target_properties( ${GlobalSearch} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}" )
    endif( MPI_LINK_FLAGS )
	
    target_link_libraries( ${GlobalSearch} ${MPI_LIBRARIES} )
	message("Set MPI FINISH")
else( ${USE_MPI} )
    message( STATUS "-- ${GlobalSearch} - NOT BUILD!"  )
endif( USE_MPI )

target_link_libraries(${GlobalSearch} ${pugixml} ${EXAMIN_MPFR_LIB} ${EXAMIN_MPIR_LIB} ${DISLIN_LIB} ${Drawing})


#add_library(${MyLib1} STATIC ${srcs} ${hdrs}) #������� ����������� ���������� � ������ �� ���������� ${MyLib1}, ����������� � ������������ CMake ����� (��������� ���������� ����������� � �����������), � �� ����������� ����� �� ���������� ${srcs} � ${hdrs}