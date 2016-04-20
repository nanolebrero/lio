######################################################################
ifeq ($(ifort),1)
FLAGS += -fpp
else
FLAGS += -cpp
endif

ifeq ($(ifort),1)
  FC   = ifort
  LIBS+= -module ../lioamber/obj
else
  FC   = gfortran
  LIBS+= -I../lioamber/obj -J../lioamber/obj
endif

LIBS += -L../g2g -lg2g
LIBS += -L../lioamber -llio-g2g

ifeq ($(magma),1)
        DEFINE   += -Dmagma
        LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif

######################################################################
liosolo: main.f
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

clean:
	rm -rf liosolo
######################################################################
