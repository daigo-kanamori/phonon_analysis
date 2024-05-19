FORTRAN := gfortran
VERSION := -std=f2003
OPTION := --free-line-length-none
OPTIMISE := -O3   -fopenmp
MLD_DIR := 
#DEBUG = -Wall -g -fbounds-check
DEBUG = 
LIB_DIR := /home/xbtss/local/include
LAPACK := -lblas -llapack
FFTW3 := -lfftw3
MDALGORITHM := BEEMAN_MODULE
MLD = $(join $(MLD_DIR),/MLD_2)
SHAREDIR := $(join $(MLD) ,/LATTICE )

################################################################################################################################################################################################################################
################################################################################################################################################################################################################################

md : MD_MODULE.o 
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG)  $(OPTIMISE) -o md.x ${ARG} *.o  -L$(LIB_DIR) $(LAPACK) $(FFTW3)

################################################################################################################################################################################################################################


################################################################################################################################################################################################################################

cond : 
	cp $(MLD)/LATTICE/MDCONDITION_MODULE.f90 . 

################################################################################################################################################################################################################################
SYSTEM_MODULE.o : 
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o SYSTEM_MODULE.o $(MLD)/SYSTEM/SYSTEM_MODULE.f90

################################################################################################################################################################################################################################

KAHAN_MODULE.o : SYSTEM_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o KAHAN_MODULE.o $(MLD)/NUMERCIAL/KAHAN_MODULE.f90

ZHEEV_MODULE.o : SYSTEM_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o ZHEEV_MODULE.o $(MLD)/NUMERCIAL/ZHEEV_MODULE.f90 

VECTOR_MODULE.o : SYSTEM_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o VECTOR_MODULE.o $(MLD)/NUMERCIAL/VECTOR_MODULE.f90

FITTING_MODULE.o : SYSTEM_MODULE.o VECTOR_MODULE.o 
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o FITTING_MODULE.o   $(MLD)/NUMERCIAL/FITTING_MODULE.f90 

################################################################################################################

IOTK_MODULE.o : SYSTEM_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o IOTK_MODULE.o   $(MLD)/STREAM/IOTK_MODULE.f90 

################################################################################################################

FFTW_MODULE.o : SYSTEM_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o FFTW_MODULE.o   $(MLD)/FFTW/FFTW_MODULE.f90 -I$(INC_DIR) 


################################################################################################################
ENSEMBLE_MODULE.o : SYSTEM_MODULE.o VECTOR_MODULE.o KAHAN_MODULE.o IOTK_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o ENSEMBLE_MODULE.o   $(MLD)/PHYSICS/ENSEMBLE_MODULE.f90 

ENSEMBLE_DATA_MODULE.o :SYSTEM_MODULE.o VECTOR_MODULE.o KAHAN_MODULE.o FITTING_MODULE.o FFTW_MODULE.o ENSEMBLE_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o ENSEMBLE_DATA_MODULE.o   $(MLD)/PHYSICS/ENSEMBLE_DATA_MODULE.f90 

PARAMETER_MODULE.o :SYSTEM_MODULE.o 
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o PARAMETER_MODULE.o   $(MLD)/PHYSICS/PARAMETER_MODULE.f90 

################################################################################################################
IBRAV_MODULE.o : SYSTEM_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o IBRAV_MODULE.o   $(MLD)/LATTICE/IBRAV_MODULE.f90 

UNITCELL_MODULE.o :SYSTEM_MODULE.o IOTK_MODULE.o VECTOR_MODULE.o IBRAV_MODULE.o KAHAN_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o UNITCELL_MODULE.o   $(MLD)/LATTICE/UNITCELL_MODULE.f90 


LATTICE_MODULE.o :SYSTEM_MODULE.o IOTK_MODULE.o PARAMETER_MODULE.o VECTOR_MODULE.o UNITCELL_MODULE.o ENSEMBLE_MODULE.o KAHAN_MODULE.o MDALGORITHM_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o LATTICE_MODULE.o   $(MLD)/LATTICE/LATTICE_MODULE.f90 

LD_MODULE.o :SYSTEM_MODULE.o IOTK_MODULE.o PARAMETER_MODULE.o UNITCELL_MODULE.o ZHEEV_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o LD_MODULE.o   $(MLD)/LATTICE/LD_MODULE.f90 

MD_MODULE.o :SYSTEM_MODULE.o PARAMETER_MODULE.o IOTK_MODULE.o LATTICE_MODULE.o ENSEMBLE_MODULE.o ENSEMBLE_DATA_MODULE.o KAHAN_MODULE.o MDCONDITION_MODULE.o MDALGORITHM_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o MD_MODULE.o   $(MLD)/LATTICE/MD_MODULE.f90  

MDALGORITHM_MODULE.o :SYSTEM_MODULE.o VECTOR_MODULE.o KAHAN_MODULE.o ENSEMBLE_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o $(MDALGORITHM).o   $(MLD)/MDALGORITHM/$(MDALGORITHM).f90 

MDCONDITION_MODULE.o : LATTICE_MODULE.o
	$(FORTRAN) $(VERSION) $(OPTION) $(DEBUG) -c $(OPTIMISE) -o MDCONDITION_MODULE.o   $(SHAREDIR)/MDCONDITION_MODULE.f90

################################################################################################################################################################################################################################
 


.PHONY: clean
clean : 
	rm *.o
	rm *.mod

