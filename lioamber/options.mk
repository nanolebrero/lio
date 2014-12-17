######################################################################
# OPTIONS
######################################################################
DEFINE   += -Dpack -fPIC

ifeq ($(print), 1)
  DEFINE += -DPRINT_MATRICES
endif

ifeq ($(profile),1)
  PROFILE  = -p
else
  PROFILE  =
endif

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

######################################################################
