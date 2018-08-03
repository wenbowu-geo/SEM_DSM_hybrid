SETUP= ../../setup/
IDIR = ../../setup/
# E : executables directory
E = ../../bin
# O : objects directory
O = ../../obj/cubedsphere_topo_forRedineOcean
# L : libraries directory
L = ../../lib
FC = ifort
FCFLAGS = #-g
FLAGS_CHECK = -O3 -assume byterecl -traceback -ftrapuv -ftz -I$(IDIR)
FLAGS_NO_CHECK = -O3 -assume byterecl -check nobounds -ftz -I$(IDIR)

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK) 
FCCOMPILE_NO_CHECK = ${FC} ${FCFLAGS} $(FLAGS_NO_CHECK) 
FCLINK = $(FCCOMPILE_NO_CHECK)

CC = gcc
CFLAGS = -g -O2 $(CPPFLAGS)
#COND_MPI_CPPFLAGS =

AR = ar
ARFLAGS = cru
RANLIB = ranlib


# objects toggled between the pure Fortran and Pyrized version
COND_PYRE_OBJECTS = $O/read_value_parameters.o $O/get_value_parameters.o

libcubedspheretopo_a_OBJECTS = \
	$O/read_parameter_file.o \
	$O/euler_angles.o \
	$O/coordi_convert.o \
	$O/cubedsphere_topo.o \
	$O/save_topo.o \
	$O/read_parameter_SEMtoTele.o \
	$O/search_topo.o \
	$(EMPTY_MACRO)


LIBCUBEDTOPO = $(COND_PYRE_OBJECTS) $L/libcubedspheretopo.a

# objects for the pure Fortran version
XCUBEDTOPO_OBJECTS = $O/program_cubedsphere_topo.o $(LIBCUBEDTOPO)

####
#### targets
####

# default targets for the pure Fortran version
DEFAULT = \
	cubedsphere_topo \
	$(EMPTY_MACRO)



default: $(DEFAULT)

all: clean default

cubedsphere_topo: xcubedsphere_topo_forRedineOcean


####
#### rules for executables
####

# rules for the pure Fortran version
xcubedsphere_topo_forRedineOcean: $(XCUBEDTOPO_OBJECTS) 
	${FCLINK} -o ${E}/xcubedsphere_topo_forRedineOcean $(XCUBEDTOPO_OBJECTS) 

# solver also depends on values from mesher



clean:
	rm -f $O/* *.o xcubedsphere_topo_forRedineOcean

###
### rule for the archive library
###

$L/libcubedspheretopo.a: $(libcubedspheretopo_a_OBJECTS)
	-rm -f $L/libcubedspheretopo.a
	$(AR) $(ARFLAGS) $L/libcubedspheretopo.a $(libcubedspheretopo_a_OBJECTS)
	$(RANLIB) $L/libcubedspheretopo.a

####
#### rule to build each .o file below
####




$O/program_cubedsphere_topo.o: program_cubedsphere_topo.f90
	${FCCOMPILE_CHECK} -c -o $O/program_cubedsphere_topo.o program_cubedsphere_topo.f90 

$O/cubedsphere_topo.o: ${SETUP}constants.h cubedsphere_topo.f90
	${FCCOMPILE_CHECK} -c -o $O/cubedsphere_topo.o cubedsphere_topo.f90 

$O/euler_angles.o: ${SETUP}constants.h euler_angles.f90
	${FCCOMPILE_CHECK} -c -o $O/euler_angles.o euler_angles.f90

$O/read_parameter_SEMtoTele.o: ${SETUP}constants.h read_parameter_SEMtoTele.f90
	${FCCOMPILE_CHECK} -c -o $O/read_parameter_SEMtoTele.o read_parameter_SEMtoTele.f90

$O/read_parameter_file.o: ${SETUP}constants.h read_parameter_file.f90
	${FCCOMPILE_CHECK} -c -o $O/read_parameter_file.o read_parameter_file.f90 

$O/read_value_parameters.o: ${SETUP}constants.h read_value_parameters.f90
	${FCCOMPILE_CHECK} -c -o $O/read_value_parameters.o read_value_parameters.f90 

$O/get_value_parameters.o: ${SETUP}constants.h get_value_parameters.f90
	${FCCOMPILE_CHECK} -c -o $O/get_value_parameters.o get_value_parameters.f90 

$O/coordi_convert.o: ${SETUP}constants.h coordi_convert.f90
	${FCCOMPILE_CHECK} -c -o $O/coordi_convert.o coordi_convert.f90 

$O/save_topo.o: ${SETUP}constants.h save_topo.f90
	${FCCOMPILE_CHECK} -c -o $O/save_topo.o save_topo.f90 

$O/search_topo.o: ${SETUP}constants.h search_topo.f90
	${FCCOMPILE_CHECK} -c -o $O/search_topo.o search_topo.f90

