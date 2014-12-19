######################################################################
#
# This file has the dependency scheme of all lio files:
#
# (1) First we list all objects to be made.
# (2) Then set the dependency for all objects that use modules.
# (3) Then set specific flags for certain objects.
#
######################################################################
# Object List
objects += liomain.o SCF.o SCFop.o TD.o
objects += dip.o dipmem.o jarz.o magnus.o predictor.o mulliken.o
objects += dft_get_mm_forces.o dft_get_qm_forces.o
objects += matmuldiag.o matmulnano.o matmulnanoc.o
objects += conmut.o conmutc.o conmutcc.o
objects += init.o init_amber.o lio_init.o lio_finalize.o
objects += alg.o drive.o func.o grid.o
objects += int1.o int1G.o int2.o int2G.o
objects += int3lu.o int3mem.o int3mems.o int3G.o
objects += intsol.o intsolG.o intsolGs.o
objects += intfld.o intSG.o
objects += FixMessRho.o get_unit.o #PackedStorage.f
objects += garcha_mod.o
objects += liokeys.o
objects += sysdata.o
objects += mathsubs.o
objects += maskrmm.o
#
######################################################################
# garcha_mod
objlist := dft_get_mm_forces.o dft_get_qm_forces.o
objlist += dip.o dipmem.o drive.o grid.o init_amber.o init.o
objlist += int1.o int2.o int3lu.o int3mem.o int3mems.o intfld.o intsol.o
objlist += int1G.o int2G.o int3G.o intSG.o intsolG.o intsolGs.o
objlist += jarz.o lio_finalize.o predictor.o SCF.o SCF_in.o SCFop.o

$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/garcha_mod.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/garcha_mod.mod
######################################################################
# Custom flags
myflags :=
ifeq ($(non_optimize),1)
  optim1 = -O0
  optim3 = -O0
else
  optim1 = -O1
  optim3 = -O3
endif

objlist := matmulnano.o matmuldiag.o conmut.o int3lu.o
$(objlist:%.o=$(obj_path)/%.o) : private myflags += $(optim3) -parallel

objlist := dip.o SCFop.o
objlist += intfld.o int1G.o int3G.o intsolG.o intsolGs.o intsol.o
$(objlist:%.o=$(obj_path)/%.o) : private myflags += $(optim1)

objlist := SCF.o TD.o ehrenfest.o magnus.o predictor.o
objlist += FixMessRho.o get_unit.o mulliken.o PackedStorage.f
objlist += conmutcc.o conmutc.o matmulnanoc.o
objlist += init_amber.o init.o lio_init.o liomain.o lio_finalize.o
objlist += dft_get_mm_forces.o dft_get_qm_forces.o
objlist += alg.o drive.o func.o grid.o dipmem.o jarz.o
objlist += int1.o int2.o int2G.o int3mem.o int3mems.o intSG.o
objlist += garcha_mod.o
$(objlist:%.o=$(obj_path)/%.o) : private myflags += $(optim3) -mp1 -i
######################################################################
