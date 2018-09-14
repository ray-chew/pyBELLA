CC = gcc

HDFLIBDIR	  = /opt/local/lib
HDFINC		  = /opt/local/include

CFLAGS		  = -Aa -I$(HDFINC) 

DEST		  = ../MakeBuilds 

EXTHDRS		  = /usr/include/math.h \
		/usr/include/stdarg.h \
		/usr/include/stdio.h \
		/usr/include/stdlib.h \
#		/usr/include/sys/stdsyms.h \
		/usr/include/sys/time.h \
		/usr/include/time.h \
#		/usr/include/varargs.h \
		/usr/include/memory.h \
		/opt/local/include/dfsd.h

HDRS		  = ../Common.h \
		../Input/userdata.h \
		../Input/boundary.h \
		../Input/enum_bdry.h \
		../Management/data.h \
		../Management/enumerator.h \
		../Management/error.h \
		../Management/io.h \
		../Management/Main.h \
		../Management/memory.h \
		../Management/variable.h \
		../Management/warning.h \
		../NumericsFundamentals/math_own.h \
		../NumericsFundamentals/Discretization/kgrid.h \
		../NumericsFundamentals/Discretization/time_discretization.h \
		../NumericsFundamentals/lapc/BiCGSTAB.h \
		../NumericsFundamentals/ofGeneralInterest/SimpleUtilities.h \
		../NumericsFundamentals/ofGeneralInterest/ThomasAlgorithmus.h \
		../Physics/Gasdynamics/Eos.h \
		../Physics/Gasdynamics/explicit.h \
		../Physics/Gasdynamics/Gasdynamics.h \
		../Physics/Gasdynamics/limiter.h \
		../Physics/Gasdynamics/numerical_flux.h \
		../Physics/Gasdynamics/recovery.h \
		../Physics/Gasdynamics/thermodynamic.h \
		../Physics/Hydrostatics/Hydrostatics.h \
		../Physics/LowMach/mpv.h \
		../Physics/LowMach/MAC-projection/flux_correction.h \
		../Physics/LowMach/MAC-projection/laplacian_cells.h \
		../Physics/LowMach/MAC-projection/variable_coefficient_poisson_cells.h \
		../Physics/LowMach/Second-projection/second_projection.h \
		../Physics/LowMach/Second-projection/laplacian_nodes.h \
		../Physics/LowMach/Second-projection/variable_coefficient_poisson_nodes.h \

HDRFLAGS	= -I../ -I../Input -I../Management -I../NumericsFundamentals \
		  -I../NumericsFundamentals/Discretization -I../NumericsFundamentals/lapc \
		  -I../NumericsFundamentals/ofGeneralInterest -I../Physics/Gasdynamics \
		  -I../Physics/Hydrostatics -I../Physics/LowMach -I../Physics/LowMach/MAC-projection \
		  -I../Physics/LowMach/Second-projection

INSTALL		= install

LD		= $(CC)

LDFLAGS		= 

LIBS		= $(HDFLIBDIR)/libmfhdf.dylib \
		$(HDFLIBDIR)/libdf.dylib  $(HDFLIBDIR)/libjpeg.a $(HDFLIBDIR)/libz.dylib -lm

MAKEFILE	= Makefile

OBJS		= userdata_InternalWave_LongWave.o \
		boundary.o \
		data.o \
		io.o \
		Main.o \
		memory.o \
		variable.o \
		kgrid.o \
		time_discretization.o \
		BiCGSTAB.o \
		SimpleUtilities.o \
		Thomas_Algorithmus_2.o \
		Eos.o \
		explicit.o \
		Gasdynamics.o \
		limiter.o \
		numerical_flux.o \
		recovery.o \
		thermodynamic.o \
		Hydrostatics.o \
		mpv.o \
		flux_correction.o \
		laplacian_cells.o \
		variable_coefficient_poisson_cells.o \
		second_projection.o \
		laplacian_nodes.o \
		variable_coefficient_poisson_nodes.o \

PRINT		  = pr

PROGRAM	   = rklm

SHELL		  = /bin/sh

SRCS		  = ../Input/userdata_InternalWave_LongWave.c \
		../Input/boundary.c \
		../Management/data.c \
		../Management/io.c \
		../Management/Main.c \
		../Management/memory.c \
		../Management/variable.c \
		../NumericsFundamentals/Discretization/kgrid.c \
		../NumericsFundamentals/Discretization/time_discretization.c \
		../NumericsFundamentals/lapc/BiCGSTAB.c \
		../NumericsFundamentals/ofGeneralInterest/SimpleUtilities.c \
		../NumericsFundamentals/ofGeneralInterest/Thomas_Algorithmus_2.c \
		../Physics/Gasdynamics/Eos.c \
		../Physics/Gasdynamics/explicit.c \
		../Physics/Gasdynamics/Gasdynamics.c \
		../Physics/Gasdynamics/limiter.c \
		../Physics/Gasdynamics/numerical_flux.c \
		../Physics/Gasdynamics/recovery.c \
		../Physics/Gasdynamics/thermodynamic.c \
		../Physics/Hydrostatics/Hydrostatics.c \
		../Physics/LowMach/mpv.c \
		../Physics/LowMach/MAC-projection/flux_correction.c \
		../Physics/LowMach/MAC-projection/laplacian_cells.c \
		../Physics/LowMach/MAC-projection/variable_coefficient_poisson_cells.c \
		../Physics/LowMach/Second-projection/second_projection.c \
		../Physics/LowMach/Second-projection/laplacian_nodes.c \
		../Physics/LowMach/Second-projection/variable_coefficient_poisson_nodes.c \

all:		$(PROGRAM)

$(PROGRAM):	 $(OBJS) 
		@echo Linking $(PROGRAM) ...
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM) $(LIBS) 
		@echo done

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:		   $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@$(MAKE) -f $(MAKEFILE) ROOT=$(ROOT) DEST=$(DEST) install

###
userdata_InternalWave_LongWave.o : 
	$(CC) -c $(HDRFLAGS) ../Input/userdata_InternalWave_LongWave.c

boundary.o : 
	$(CC) -c $(HDRFLAGS) ../Input/boundary.c

data.o : 
	$(CC) -c $(HDRFLAGS) ../Management/data.c

io.o : 
	$(CC) -c $(HDRFLAGS) ../Management/io.c

Main.o : 
	$(CC) -c $(HDRFLAGS) ../Management/Main.c

memory.o : 
	$(CC) -c $(HDRFLAGS) ../Management/memory.c

variable.o : 
	$(CC) -c $(HDRFLAGS) ../Management/variable.c

kgrid.o : 
	$(CC) -c $(HDRFLAGS) ../NumericsFundamentals/Discretization/kgrid.c

time_discretization.o : 
	$(CC) -c $(HDRFLAGS) ../NumericsFundamentals/Discretization/time_discretization.c

BiCGSTAB.o : 
	$(CC) -c $(HDRFLAGS) ../NumericsFundamentals/lapc/BiCGSTAB.c

SimpleUtilities.o : 
	$(CC) -c $(HDRFLAGS) ../NumericsFundamentals/ofGeneralInterest/SimpleUtilities.c

Thomas_Algorithmus_2.o : 
	$(CC) -c $(HDRFLAGS) ../NumericsFundamentals/ofGeneralInterest/Thomas_Algorithmus_2.c

Eos.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/Gasdynamics/Eos.c

explicit.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/Gasdynamics/explicit.c

Gasdynamics.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/Gasdynamics/Gasdynamics.c

limiter.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/Gasdynamics/limiter.c

numerical_flux.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/Gasdynamics/numerical_flux.c

recovery.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/Gasdynamics/recovery.c

thermodynamic.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/Gasdynamics/thermodynamic.c

Hydrostatics.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/Hydrostatics/Hydrostatics.c

mpv.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/LowMach/mpv.c

flux_correction.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/LowMach/MAC-projection/flux_correction.c
	
laplacian_cells.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/LowMach/MAC-projection/laplacian_cells.c

variable_coefficient_poisson_cells.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/LowMach/MAC-projection/variable_coefficient_poisson_cells.c

second_projection.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/LowMach/Second-projection/second_projection.c

laplacian_nodes.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/LowMach/Second-projection/laplacian_nodes.c

variable_coefficient_poisson_nodes.o : 
	$(CC) -c $(HDRFLAGS) ../Physics/LowMach/Second-projection/variable_coefficient_poisson_nodes.c

