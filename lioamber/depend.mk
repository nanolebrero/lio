######################################################################
# FILE DEPENDENCIES
######################################################################
# OBJECT LIST
objects := $(objects) alg.o drive.o func.o grid.o
objects := $(objects) conmutcc.o conmutc.o conmut.o
objects := $(objects) dft_get_mm_forces.o dft_get_qm_forces.o
objects := $(objects) dip.o dipmem.o
objects := $(objects) ehrenfest.o
objects := $(objects) FixMessRho.o
objects := $(objects) get_unit.o
objects := $(objects) init_amber.o
objects := $(objects) init.o
objects := $(objects) int1.o int1G.o
objects := $(objects) int2.o int2G.o
objects := $(objects) int3mem.o int3mems.o int3lu.o int3G.o
objects := $(objects) intfld.o
objects := $(objects) intSG.o
objects := $(objects) intsol.o intsolG.o intsolGs.o
objects := $(objects) jarz.o
objects := $(objects) lio_finalize.o
objects := $(objects) lio_init.o
objects := $(objects) liomain.o
objects := $(objects) magnus.o predictor.o
objects := $(objects) matmuldiag.o matmulnanoc.o matmulnano.o
objects := $(objects) mulliken.o
objects := $(objects) SCF.o SCFop.o
objects := $(objects) TD.o

# GARCHAMOD DEPENDS
objlist :=
objlist += dft_get_mm_forces.o dft_get_qm_forces.o
objlist += dip.o dipmem.o drive.o grid.o init_amber.o init.o
objlist += int1.o int2.o int3lu.o int3mem.o int3mems.o intfld.o intsol.o
objlist += int1G.o int2G.o int3G.o intSG.o intsolG.o intsolGs.o
objlist += jarz.o lio_finalize.o predictor.o SCF.o SCF_in.o SCFop.o
$(objlist) : garcha_mod.o

######################################################################
