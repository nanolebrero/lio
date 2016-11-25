######################################################################
# MAKE DEPEND
#
# This file contains the information of dependencies between the
# different objects that make up LIO. The list "objects" contains
# ALL the objects to be linked for making LIO, including both free
# subroutines and modules (their corresponding .o files).
#
# Each module section should:
# 
#  (*) Include in the list of objects the object file that contains
#      the module and all special objects needed.
#
#  (*) The path where to find all sources to compile the objects
#      that were included in the list. This means the .f sources
#      and the .mk that clarify internal dependencies.
#
#  (*) Include in this makefile any file where the internal
#      dependencies of the module are clarified.
#
#  (*) Specify the dependencies related to this object (it is
#      always preferable to specify which other files depend
#      on that module that the other way around).
#
# Please note that the obj_path is explicitly set in each case; 
# vpath can't be trusted with handling files generated during
# compiling.
#
######################################################################
objects += liomain.o SCF.o SCFop.o SCF_in.o TD.o cubegen.o
objects += dip.o dipmem.o jarz.o magnus.o predictor.o mulliken.o
objects += dft_get_mm_forces.o dft_get_qm_forces.o
objects += matmuldiag.o fock_commuts.o
objects += init_lio.o lio_finalize.o
objects += drive.o func.o grid.o
objects += int1.o int1G.o int2.o int2G.o
objects += int3lu.o int3mem.o  int3G.o
objects += intsol.o intsolG.o intsolGs.o
objects += intfld.o intSG.o
objects += FixMessRho.o PackedStorage.o
objects += liokeys.o
objects += sysdata.o
objects += ECP_mod.o
objects += test_subs.o
objects += esp_funct.o
objects += readECP.o
objects += intECP.o
objects += density.o
objects += extras.o
objects += fterm_biaspot.o lowdinpop.o
objects += elec.o
objects += SCF_gro.o
objects += transport.o
# objects += init_gromacs.o init_amber.o init.o lio_init.o
#
#
#     Trying a new way of makefile organization: put every important
# information inside of the module.mk
######################################################################
# linear_algebra
objects += linear_algebra.o
src_paths += linear_algebra
include linear_algebra/linear_algebra.mk

# liosubs
objects += liosubs.o
src_paths += liosubs
include liosubs/liosubs.mk



# garcha_mod: Description pending
######################################################################
objects   += garcha_mod.o
src_paths += liomods
include liomods/liomods.mk

tmplist := dft_get_mm_forces.o dft_get_qm_forces.o
tmplist += dip.o dipmem.o drive.o grid.o init_amber.o init.o
tmplist += int1.o int2.o int3lu.o int3mem.o intfld.o intsol.o
tmplist += int1G.o int2G.o int3G.o intSG.o intsolG.o intsolGs.o
tmplist += jarz.o lio_finalize.o predictor.o
tmplist += SCF.o SCF_in.o SCFop.o TD.o cubegen.o
tmplist += maskrmm.o
$(tmplist:%.o=$(obj_path)/%.o) : $(obj_path)/garcha_mod.mod


# maskrmm: Description pending
######################################################################
objects   += maskrmm.o
src_paths += maskrmm
include maskrmm/maskrmm.mk


# mathsubs: Description pending
######################################################################
objects   += mathsubs.o
src_paths += mathsubs
include mathsubs/mathsubs.mk

tmplist := SCF.o SCFop.o
$(tmplist:%.o=$(obj_path)/%.o) : $(obj_path)/mathsubs.mod


# general_module: Description pending
######################################################################
objects   += general_module.o
src_paths += general_module
include general_module/general_module.mk

tmplist := SCF.o
$(tmplist:%.o=$(obj_path)/%.o) : $(obj_path)/general_module.mod


# ECP_mod
objlist := SCF.o
objlist += readECP.o
objlist += drive.o
objlist += intECP.o
objlist += init.o extras.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/ECP_mod.mod

# intECP
objlist := SCF.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/intECP.o

objlist := drive.o
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/readECP.o

#esp_funct
objlist := intECP.o
objlist += ECP_mod.mod
$(objlist:%.o=$(obj_path)/%.o) : $(obj_path)/esp_funct.mod

# cublas: Description pending
######################################################################
ifeq ($(cublas),1)
objects   += cublasmath.o fortran.o
src_paths += cublasmath
include cublasmath/cublasmath.mk

tmplist := cublasmath.o
$(tmplist:%.o=$(obj_path)/%.o) : $(obj_path)/garcha_mod.mod
$(tmplist:%.o=$(obj_path)/%.o) : $(obj_path)/fortran.o

tmplist := SCF.o SCFop.o TD.o
$(tmplist:%.o=$(obj_path)/%.o) : $(obj_path)/cublasmath.mod

$(obj_path)/fortran.o: $(CUDA_HOME)/src/fortran.c
	$(CC) -fPIC -DCUBLAS_GFORTRAN -O3 -c $< -o $@ -I$(CUDA_HOME)/include
endif


######################################################################
