#Definition of the compiler
CC=gfortran
#Specific flags
# -pg is for gprofiling which tells you how much time you spend on which function
# -no-pie is also part of gprofiling means position independent and tells to compiler not to generate the pie version of the code
# -fno-builtin is part of gprofiling whcih tells not to try optimise my code by replacing my functions with builtin functions or inline functions
# -O is for debugging. higher number after O gives less debugging layers
FLAGS= -g -O0 -Wno-implicit -pg -no-pie -fno-builtin -freal-8-real-16 #-Wall
#Sources
SOURCES= ./precision.F90 ./Structures.F90 ./Generic.F90 ./strings.F90 ./evaluate.F90 ./Msh2Tri.F90 ./structured_meshgen.F90 ./ShapFun.F90 ./ShapFun_unstruc.F90 ./splitting.F90 ./matrices.F90  ./get_vtk_files.F90 ./transport_rect.F90 ./transport_tri.F90 ./transport_tri_unstr.F90 ./transport_tri_semi.F90 ./amin.F90 ./main.F90
#Souces compiled
CSOURCES = ./precision.o ./Structures.o ./Generic.o ./strings.o ./evaluate.o ./Msh2Tri.o ./structured_meshgen.o ./ShapFun.o ./ShapFun_unstruc.o ./splitting.o ./matrices.o ./get_vtk_files.o ./transport_rect.o ./transport_tri.o ./transport_tri_unstr.o ./transport_tri_semi.o ./amin.o ./main.o
#Dependencies
libraries= #-Wl,-rpath,/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/lib -Wl,-rpath,/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/lib -Wl,-rpath,/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/lib -lpetsc
#Include files
include = #-I/usr/lib/petscdir/3.8.3/include -I/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/include
#myprog:
#Name of the executable
EXECUTABLE=runmycase



all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(SOURCES) 
	@echo 'Building target: $@'
	@echo 'Invoking: Fortran Linker'
	$(CC) -c $(FLAGS) $(include) $(SOURCES) 
	$(CC) -o $(EXECUTABLE) $(FLAGS) $(CSOURCES) $(libraries) 
	@echo 'Finished building target: $@'
	@echo ' '

	#rm -r *.o *.mod


# just change the name of the files for SCR
# to run this first type make
# then, $ ./tri_DG
# compile
#FC = gfortran
#CFLAGS = -std=f2008 -O -Wall -fcheck=all -g -fbacktrace

#testfile : Tri_stab2.f90
#	$(FC) $(CFLAGS) -o testrun Tri_stab2.f90

#clean:
#	rm -r testfile *.o *.mod


# https://www.youtube.com/watch?v=OMiGbKS00mU&list=PLOU8LxhyFylLS298Sea2-gYvO5Lj8HZsP&index=2
# # cmpiler flags
# #standards
# # CFLAGS = -std=f008ts
# # warning flag
# CFLAGS += -Wall
# #debgging options
# CFLAGS += -fPIC -fmax-errors=3
#
# #source file
# SRCS = test.f90
# OBJs = $(SRCS:=.o)
#
# # excecutable
# MAIN = main
#
# # cmpile Project
# all : (MAIN)
# 	@echo Model compiled
#
# $(MAIN) : $(OBJS)
# 	$(FC) $(CFLAGS) -o $(MAIN) $(OBJS)
#
# .SUFFIXES : .o .f90
#
# .f90.o :
# 	$(FC) $(CFLAGS) -c $<
#
# clean :
# 	$(RM) *.o *.mod $(MAIN)
