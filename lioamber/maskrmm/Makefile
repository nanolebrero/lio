######################################################################
# COMPILER AND FLAGS
FC     = ifort
FFLAGS = -g -O3 -c -cpp -Dpack -fPIC
#
#
# DEPENDENCIES
INCLUDES  =
INCLUDES += spunpack.f sprepack.f
INCLUDES += rhomess.f  rhomess_h.f
INCLUDES += rhofix.f   rhofix_h.f
#
#
# MODULE COMPILATION COMMANDS
all : maskrmm.o
	cp maskrmm.mod ../.
	cp maskrmm.o ../obj/.


maskrmm.o : $(INCLUDES)

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@
#
#
# CLEAN
.PHONY: clean
clean :
	rm -f maskrmm.o maskrmm.mod
######################################################################
