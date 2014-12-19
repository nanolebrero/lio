######################################################################
# LIBLIO MAKEFILE
######################################################################
src_paths := liomods
src_paths += maskrmm
src_paths += mathsubs
obj_path  := obj
vpath %.f   $(src_paths)
vpath %.mod $(obj_path)
vpath %.o   $(obj_path)

include make_options.mk

all: liblio-g2g.so

objects :=
include make_depend.mk
include maskrmm/maskrmm.mk
include mathsubs/mathsubs.mk

liblio-g2g.so: $(objects:%.o=$(obj_path)/%.o)
	$(FC) -shared -fPIC -o $@ $^ $(LIBS)

include make_rules.mk
$(objects:%.o=$(obj_path)/%.o) : Makefile
$(objects:%.o=$(obj_path)/%.o) : make_rules.mk
$(objects:%.o=$(obj_path)/%.o) : make_depend.mk
$(objects:%.o=$(obj_path)/%.o) : make_options.mk
######################################################################
