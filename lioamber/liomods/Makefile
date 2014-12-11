#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# SOURCES AND PATHS
#--------------------------------------------------------------------#
MODPATH = ../mod
OBJPATH = ../obj
SRCS    = $(wildcard *.f)
OBJS    = $(SRCS:%.f=%.o)
MODS    = $(SRCS:%.f=%.mod)
#
#
# COMPILING ORDERS
#--------------------------------------------------------------------#
FC      = ifort
FFLAGS  = -g -O3 -c -cpp -Dpack -fPIC -module $(MODPATH)

all : objects modules
objects : $(OBJS:%.o=$(OBJPATH)/%.o)
modules : $(MODS:%.mod=$(MODPATH)/%.mod)

$(OBJPATH)/%.o   : %.f ../obj
	$(FC) $(FFLAGS) -c $< -o $@
../obj:
	mkdir ../obj

$(MODPATH)/%.mod : %.f ../mod
	$(FC) $(FFLAGS) -c $< -o $(@:$(MODPATH)/%.mod=$(OBJPATH)/%.o)
../mod:
	mkdir ../mod
#
#
# CLEAN
#--------------------------------------------------------------------#
.PHONY: clean
clean :
	rm -f $(OBJS:%.o=$(OBJPATH)/%.o) $(MODS:%.mod=$(MODPATH)/%.mod)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
