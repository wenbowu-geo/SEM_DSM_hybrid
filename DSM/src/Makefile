FC      = mpif77
FC90    = ifort
CC	= mpicc
FFLAGS  = -O 
#-check nobounds
CFLAGS	= $(FFLAGS) 
OBJS    = param.o matrix.o interpolate.o compwave.o msub.o savespec.o check_add_node.o dsmti.o



all:
	make dsmti

dsmti: ${OBJS}
	$(FC) $(OBJS) $(FFLAGS) -o dsmti

%.o: %.f90
	$(FC90) -c $(FFLAGS) $<
%.o: %.f
	$(FC) -c $(FFLAGS) $<
.c.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f $(OBJS) dsmti *~
