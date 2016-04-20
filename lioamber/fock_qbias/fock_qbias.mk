######################################################################
# EXTERNAL DEPENDENCIES : Which other objets use this module? .o

CLIENTS := SCF.o TD.o
$(CLIENTS:%.o=$(obj_path)/%.o) : $(obj_path)/fock_qbias.o

######################################################################
# INTERNAL DEPENDENCIES : Which files does this module use? .f

SOURCES := fock_qbias.mk
SOURCES += fock_qbias.f90
SOURCES += fqbias_setup.f90
SOURCES += read_atom_prop.f90
SOURCES += align_ab_prop.f90
SOURCES += fqbias_calc.f90
SOURCES += gaussian_shape.f90
SOURCES += fqbias_lowdinpop.f90
$(obj_path)/fock_qbias.o : $(SOURCES)

######################################################################
