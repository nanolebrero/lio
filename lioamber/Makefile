######################################################################
# LIBLIO MAKEFILE
######################################################################
#
#
#
#
######################################################################
all: liblio-g2g.so

mod_path  := mod
obj_path  := obj
src_paths := liomods maskrmm mathsubs
vpath %.mod $(mod_path)
vpath %.o   $(obj_path)
vpath %.f   $(src_paths)

# Flags
FC     = ifort
FFLAGS =
FFLAGS += -g -module $(mod_path)
FFLAGS += $(DEFINE) $(PROFILE)
DEFINE += -Dpack -fPIC




# Libraries
LIBS_PATH =
LIBS_PATH += -L../g2g -L/usr/lib -L/usr/lib64
LIBS_LIBS += -lg2g -lstdc++


# MKL libs and mods
LIBS_PATH += -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include
LIBS_LIBS += -lmkl_lapack95_lp64
LIBS_LIBS += -lmkl_intel_lp64
LIBS_LIBS += -lmkl_intel_thread
LIBS_LIBS += -lmkl_core


LIBS_MOAR = -liomp5 -lpthread -lm -Wl,-rpath='$$ORIGIN/' -Wl,-rpath='$$ORIGIN/../g2g'

LIBS = $(LIBS_PATH) $(PROFILE) $(LIBS_LIBS) $(LIBS_MOAR)
#
#
#
######################################################################

CFLAGS   := -mp1 -ip -cpp  $(DEFINE) $(PROFILE)
CFLAGS2  := $(OPTIMIZE) -cpp  $(DEFINE) $(PROFILE)
CFLAGS3  := $(OPTIMIZE) -cpp -parallel   $(DEFINE) $(PROFILE)

ifeq ($(non_optimize),1)
  CFLAGS  += -O0
  CFLAGS2 += -O0
  CFLAGS3 += -O0
else
  CFLAGS  += -O3
  CFLAGS2 += -O1
  CFLAGS3 += -O3
endif

ifeq ($(magma),1)
	DEFINE   += -Dmagma
	LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif


## Targets: make ##

all: dirs liblio-g2g.so

dirs:
	mkdir -p obj
	mkdir -p mod

liblio-g2g.so: dirs $(GARCHA_OBJ:%.o=obj/%.o)
	ifort -shared -fPIC -o $@ $(GARCHA_OBJ:%.o=obj/%.o) $(LIBS)


## Targets: clean ##
clean:
	rm -rf obj
	rm -rf libgarcha-g2g.so liblio-g2g.so


# Different flags for some .o
obj/intfld.o: intfld.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/dipG.o: dipG.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/dip.o: dip.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/int1G.o: int1G.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/int3G.o: int3G.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/intsolG.o: intsolG.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS)  $<

obj/init.o: init.f
	ifort -o $@ $(CFLAGS) $(G2G_CFLAGS)  $<

obj/SCFop.o: SCFop.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/intsolGs.o: intsolGs.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS) $<

obj/matmulnano.o: matmulnano.f
	ifort -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<

obj/conmut.o: conmut.f
	ifort -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<

obj/matmuldiag.o: matmuldiag.f
	ifort -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<

obj/int3lu.o: int3lu.f
	ifort -o $@ $(CFLAGS3)  $(G2G_CFLAGS) $<

obj/intsol.o: intsol.f
	ifort -o $@ $(CFLAGS2) $(G2G_CFLAGS)  $<

# Directories dependencies
%.o   : | $(obj_path)
%.mod : | $(mod_path)
$(obj_path) $(mod_path) :
	mkdir -p $@

# Objects
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $(obj_path)/$@

.PHONY: clean
clean:
	rm -rf obj mod liblio-g2g.so

######################################################################
