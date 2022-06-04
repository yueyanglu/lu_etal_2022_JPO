#
# Makefile for HYCOM post-processing synthetic floats
#
# --- usage: make lagrangian ARCH=sun >& Make_lagrangian
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
	@exit 1
.f:
	@echo "Must have an explicit rule for" $*
	@exit 1
.c:
	@echo "Must have an explicit rule for" $*
	@exit 1

#include ../../config/$(ARCH)_setup
include $(ARCH)_setup

#
# ---------------------------------------------------------------------
# Standard part of the Makefile starts here
# ---------------------------------------------------------------------
#

all:
	/bin/csh Make_all.com

clean:
	/bin/rm *.o *.a *.mod M*log *tspt *port *port2 *vel *mn *mn2 *2p0

MODS =	mod_tracer.o mod_za.o mod_xc.o wtime.o
OBJS =	blkin.o geopar.o getdat.o zh.o matinv.o 
#fct3d.o

#
# --- executables.
#

tracers:      $(MODS) tracers.o            $(OBJS)
	$(LD) $(LDFLAGS) tracers.o    $(MODS) $(OBJS) -o tracers

#
# --- explicit dependencies for each .[fF] file.
#

geopar.o:        geopar.f      $(MODS)
getdat.o:        getdat.f      $(MODS)
#advem.o:         advem.f
#fct3d.o:         fct3d.f
matinv.o:        matinv.f       
blkin.o:         blkin.f
wtime.o:         wtime.F

mod_tracer.o: mod_tracer.F mod_xc.o mod_za.o
mod_xc.o:     mod_xc.F   
mod_za.o:     mod_za.F    mod_xc.o

zh.o:         zh.F zh_sun.f zh_t3e.f zh_ymp.f
