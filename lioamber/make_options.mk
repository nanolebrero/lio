######################################################################
# Compiler and common flags
FC     = ifort
FFLAGS =
FFLAGS += -Dpack -fPIC
FFLAGS += -g -fpp -module $(obj_path)
FFLAGS += $(myflags) $(DEFINE) $(PROFILE)

# Libs and mods (common and MKL)
LIBS += -L../g2g -L/usr/lib -L/usr/lib64
LIBS += -lg2g -lstdc++
LIBS += -liomp5 -lpthread -lm -Wl,-rpath='$$ORIGIN/' -Wl,-rpath='$$ORIGIN/../g2g'
LIBS += -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include
LIBS += -lmkl_lapack95_lp64 -lmkl_intel_lp64
LIBS += -lmkl_intel_thread -lmkl_core

ifeq ($(print), 1)
  DEFINE += -DPRINT_MATRICES
endif

ifeq ($(profile),1)
  PROFILE = -p
else
  PROFILE =
endif

ifeq ($(magma),1)
  DEFINE += -Dmagma
  LIBS   += -L$(MAGMAROOT)/lib -lmagma
endif

######################################################################
