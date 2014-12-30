######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += spunpack.f sprepack.f
INCLUDES += rhomess.f  rhomess_h.f
INCLUDES += rhofix.f   rhofix_h.f
INCLUDES += rmmget.f  rmmput.f

$(obj_path)/maskrmm.o : $(INCLUDES) maskrmm.mk
######################################################################
