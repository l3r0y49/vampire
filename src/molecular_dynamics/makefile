#--------------------------------------------------------------
#          Makefile for molecular_dynamics module
#--------------------------------------------------------------

# List module object filenames
molecular_dynamics_objects =\
data.o \
initialize.o \
interface.o \
compute.o \
gen_coords.o \
neighbour_list_manage.o \
evolve_md.o \
md_main.o \
terminate.o


# Append module objects to global tree
OBJECTS+=$(addprefix obj/molecular_dynamics/,$(molecular_dynamics_objects))
