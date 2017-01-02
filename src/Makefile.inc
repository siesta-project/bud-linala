LA_src_DIR = $(TOP_DIR)/src
# Append to VPATH
VPATH := $(VPATH):$(LA_src_DIR)

BLAS_OBJECTS = \
	bud_BLAS1.o \
	bud_BLAS2.o \
	bud_BLAS3.o \
	bud_BLAS.o
$(BLAS_OBJECTS:.o=.f90): BLAS.inc
LA_src_OBJECTS += $(BLAS_OBJECTS)
bud_BLAS1.f90: BLAS1.F90 BLAS1.inc
	$(FPP) $(FPPFLAGS) $(INCLUDES) $< > $@
bud_BLAS2.f90: BLAS2.F90 BLAS2.inc
	$(FPP) $(FPPFLAGS) $(INCLUDES) $< > $@
bud_BLAS3.f90: BLAS3.F90 BLAS3.inc
	$(FPP) $(FPPFLAGS) $(INCLUDES) $< > $@
bud_BLAS.f90: BLAS.F90
	$(FPP) $(FPPFLAGS) $(INCLUDES) $< > $@


# Define objects that needs compilation
LA_src_SOURCES = $(patsubst %.o, %.f90, $(LA_src_OBJECTS))
# FPP_SOURCES are those sources that needs deletion upon make clean
LA_src_FPP_SOURCES = $(LA_src_SOURCES)
LA_src_MODS = $(patsubst %.f90, %.mod, $(LA_src_SOURCES))

SOURCES += $(LA_src_SOURCES)
FPP_SOURCES += $(LA_src_FPP_SOURCES)
F_MODS  += $(LA_src_MODS)
OBJECTS += $(LA_src_OBJECTS)


# Define rules for all objects
$(LA_src_FPP_SOURCES): INCLUDES += -I$(LA_src_DIR)


source-la-src: $(LA_src_SOURCES)


clean-la-src:
	-rm -f $(LA_src_FPP_SOURCES)
	-rm -f $(LA_src_OBJECTS)
	-rm -f *.mod


clean: clean-la-src


# Local Variables:
#  mode: makefile-gmake
# End:
