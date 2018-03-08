
#--------------------------------
# insert name of your executable
#--------------------------------
PROGRAM = mandelseq

#-------------------------------------------------
# insert compiler: 
#        icc         (serial)
#        mpicc       (mpi)
#        icc -openmp (openmp)
#-------------------------------------------------
CC  = mpicc


OBJ     = ${PROGRAM:=.o}

ARCH = _LINUX
INCLUDE_DIR = -I../Tools
LIBS =  -L../Tools -lppmwrite -lm
USER_CFLAGS =  $(INCLUDE_DIR) -D${ARCH} -g 

USER_FFLAGS =  $(INCLUDE_DIR) -r8  -g 
srcdir=.

VPATH=.:$(srcdir)

$(PROGRAM): $(OBJ) 
	$(CC) $(USER_CFLAGS) -o $(PROGRAM) $(OBJ) $(LIBS)

.c.o:
	$(CC) $(USER_CFLAGS) -c $<

clean: 
	rm -f *.o $(PROGRAM)


