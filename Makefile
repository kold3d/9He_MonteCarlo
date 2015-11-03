CDIR    =       $(CERN)
CLIBDIR =       $(CDIR)/$(CERN_LEVEL)/lib
ULIB    =       /usr/lib
XULIB    =      /usr/X11R6/lib
CFLAG	=	-c
COPT	=	-w
OPTIM	=	-O3
LINKER  =       gfortran
COMP    =       gfortran
DEBUG   =       -g
KUIPC   =       kuipc
DBGOPT  =       -Wall


LIBFLS  =  -lmathlib -lgeant -lgeant321 -lgraflib -lX11 -lgrafX11\
           -lkernlib -lpacklib -lc -lgeant321 -lgeant\
           -lphtools
            

OBJ     =  main.o uginit.o gustep.o ugeom.o uglast.o gutrev.o guout.o\
           gukine.o mykin.o getweight.o

all: mgeant

mgeant: $(OBJ)
	$(LINKER) -o mgeant $(OBJ) \
	-L$(CLIBDIR) -L$(XULIB) -L$(ULIB) $(LIBFLS)

.SUFFIXES: .for

.for.o:
	$(COMP) $(OPTIM) -c -fno-automatic -o $@ $<

clean:
	rm *.o mgeant
